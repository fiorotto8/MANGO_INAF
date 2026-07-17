#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kRadToDeg = 180.0/kPi;

struct Direction {
  double ux;
  double uy;
  double uz;
};

struct AngleValues {
  double theta;
  double phi;
  double sinTheta;
  double absUx;
};

AngleValues Angles(const Direction& direction)
{
  const double ux = std::clamp(direction.ux, -1.0, 1.0);
  const double theta = std::acos(ux);
  return {theta, std::atan2(direction.uy, direction.uz),
          std::sin(theta), std::abs(direction.ux)};
}

bool NearlyEqual(double actual, double expected, double tolerance = 1.e-12)
{
  return std::abs(actual - expected) <= tolerance;
}

int SelfTest()
{
  struct Test {
    const char* name;
    Direction direction;
    double thetaDeg;
    double phiDeg;
    bool checkPhi;
  };
  const double s3 = std::sqrt(1.0/3.0);
  const std::array<Test, 7> tests = {{
    {"+x", { 1., 0., 0.},   0.,          0., false},
    {"-x", {-1., 0., 0.}, 180.,          0., false},
    {"+y", { 0., 1., 0.},  90.,         90., true},
    {"-y", { 0.,-1., 0.},  90.,        -90., true},
    {"+z", { 0., 0., 1.},  90.,          0., true},
    {"-z", { 0., 0.,-1.},  90.,        180., true},
    {"(1,1,1)/sqrt(3)", {s3,s3,s3}, 54.7356103172453, 45., true}
  }};

  bool passed = true;
  std::cout << std::fixed << std::setprecision(12);
  for (const auto& test : tests) {
    const auto values = Angles(test.direction);
    const double thetaDeg = values.theta*kRadToDeg;
    const double phiDeg = values.phi*kRadToDeg;
    const bool rowPassed =
      NearlyEqual(thetaDeg, test.thetaDeg, 1.e-10) &&
      (!test.checkPhi || NearlyEqual(phiDeg, test.phiDeg, 1.e-10));
    passed = passed && rowPassed;
    std::cout << "ANGLE_TEST " << test.name
              << " theta_deg=" << thetaDeg
              << " phi_deg=" << phiDeg
              << " sin_theta=" << values.sinTheta
              << " abs_ux=" << values.absUx
              << " status=" << (rowPassed ? "PASS" : "FAIL") << '\n';
  }
  std::cout << "ANGLE_TEST_SUMMARY " << (passed ? "PASS" : "FAIL") << '\n';
  return passed ? 0 : 1;
}

double Quantile(std::vector<double> values, double probability)
{
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(values.begin(), values.end());
  const double position = probability*(values.size() - 1);
  const auto lower = static_cast<std::size_t>(std::floor(position));
  const auto upper = static_cast<std::size_t>(std::ceil(position));
  const double fraction = position - lower;
  return values[lower]*(1. - fraction) + values[upper]*fraction;
}

double Fraction(const std::vector<double>& values,
                const std::function<bool(double)>& predicate)
{
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  const auto count = std::count_if(values.begin(), values.end(), predicate);
  return static_cast<double>(count)/values.size();
}

struct Summary {
  long long electronCount = 0;
  long long eventCount = 0;
  double thetaMedianDeg = std::numeric_limits<double>::quiet_NaN();
  double thetaP16Deg = std::numeric_limits<double>::quiet_NaN();
  double thetaP84Deg = std::numeric_limits<double>::quiet_NaN();
  double meanSinTheta = std::numeric_limits<double>::quiet_NaN();
  double rmsSinTheta = std::numeric_limits<double>::quiet_NaN();
  double fracPlane10 = std::numeric_limits<double>::quiet_NaN();
  double fracPlane20 = std::numeric_limits<double>::quiet_NaN();
  double fracPlane30 = std::numeric_limits<double>::quiet_NaN();
  double fracSin025 = std::numeric_limits<double>::quiet_NaN();
  double fracSin050 = std::numeric_limits<double>::quiet_NaN();
  double fracSin075 = std::numeric_limits<double>::quiet_NaN();
};

Summary Summarize(const std::vector<double>& theta,
                  const std::vector<double>& sinTheta,
                  const std::set<int>& events)
{
  Summary result;
  result.electronCount = static_cast<long long>(theta.size());
  result.eventCount = static_cast<long long>(events.size());
  if (theta.empty()) return result;

  std::vector<double> thetaDeg;
  thetaDeg.reserve(theta.size());
  for (double value : theta) thetaDeg.push_back(value*kRadToDeg);
  result.thetaMedianDeg = Quantile(thetaDeg, 0.50);
  result.thetaP16Deg = Quantile(thetaDeg, 0.16);
  result.thetaP84Deg = Quantile(thetaDeg, 0.84);

  double sum = 0.;
  for (double value : sinTheta) sum += value;
  result.meanSinTheta = sum/sinTheta.size();
  double variance = 0.;
  for (double value : sinTheta) {
    const double delta = value - result.meanSinTheta;
    variance += delta*delta;
  }
  result.rmsSinTheta = std::sqrt(variance/sinTheta.size());

  result.fracPlane10 = Fraction(thetaDeg,
    [](double value) { return std::abs(value - 90.) < 10.; });
  result.fracPlane20 = Fraction(thetaDeg,
    [](double value) { return std::abs(value - 90.) < 20.; });
  result.fracPlane30 = Fraction(thetaDeg,
    [](double value) { return std::abs(value - 90.) < 30.; });
  result.fracSin025 = Fraction(sinTheta,
    [](double value) { return value < 0.25; });
  result.fracSin050 = Fraction(sinTheta,
    [](double value) { return value < 0.50; });
  result.fracSin075 = Fraction(sinTheta,
    [](double value) { return value < 0.75; });
  return result;
}

struct Record {
  int eventID = 0;
  int trackID = 0;
  int parentID = 0;
  int gasEntered = 0;
  int fullyContained = 0;
  int selected = 0;
  double generatedEnergy = 0.;
  double generatedTheta = 0.;
  double generatedSinTheta = 0.;
  double entryEnergy = 0.;
  double depositedEnergy = 0.;
  double theta = 0.;
  double phi = 0.;
  double sinTheta = 0.;
  double ux = 0.;
  double uy = 0.;
  double uz = 0.;
};

struct Accumulator {
  std::vector<double> theta;
  std::vector<double> sinTheta;
  std::set<int> events;
};

enum class Stage { Generated = 0, GasEntry = 1, Contained = 2, Selected = 3 };

const char* StageName(Stage stage)
{
  switch (stage) {
    case Stage::Generated: return "generated";
    case Stage::GasEntry: return "gas_entry";
    case Stage::Contained: return "contained";
    case Stage::Selected: return "selected";
  }
  return "unknown";
}

const char* EnergyDefinition(Stage stage)
{
  return stage == Stage::Generated ? "generated_kinetic_energy"
                                   : "deposited_energy";
}

bool Passes(Stage stage, const Record& record)
{
  switch (stage) {
    case Stage::Generated: return true;
    case Stage::GasEntry: return record.gasEntered;
    case Stage::Contained: return record.gasEntered && record.fullyContained;
    case Stage::Selected: return record.gasEntered && record.selected;
  }
  return false;
}

bool InOpenBin(double energy, double low, double high)
{
  return energy > low && energy < high;
}

void WriteNumber(std::ostream& stream, double value)
{
  if (std::isfinite(value)) stream << value;
  else stream << "NA";
}

}  // namespace

int main(int argc, char** argv)
{
  try {
    if (argc == 2 && std::string(argv[1]) == "--self-test") return SelfTest();
    if (argc != 3) {
      std::cerr << "Usage: " << argv[0]
                << " --self-test\n       " << argv[0]
                << " RAW_GEANT4_ROOT OUTPUT_DIRECTORY\n";
      return 2;
    }

    const std::filesystem::path inputPath = argv[1];
    const std::filesystem::path outputDirectory = argv[2];
    std::filesystem::create_directories(outputDirectory);
    const auto rootPath = outputDirectory/"sr90_acceptance.root";
    const auto csvPath = outputDirectory/"sr90_acceptance.csv";
    const auto pngPath = outputDirectory/"sr90_acceptance.png";

    TFile input(inputPath.c_str(), "READ");
    if (input.IsZombie()) throw std::runtime_error("cannot open input ROOT file");
    auto* tree = dynamic_cast<TTree*>(input.Get("Acceptance"));
    if (!tree) throw std::runtime_error("input ROOT file has no Acceptance tree");
    const long long inputRows = tree->GetEntries();

    Record record;
    double generatedUx = 0.;
    double generatedUy = 0.;
    double generatedUz = 0.;
    double generatedPhi = 0.;
    double generatedAbsUx = 0.;
    double absUx = 0.;
    tree->SetBranchAddress("EventID", &record.eventID);
    tree->SetBranchAddress("TrackID", &record.trackID);
    tree->SetBranchAddress("ParentID", &record.parentID);
    tree->SetBranchAddress("GeneratedEnergy_keV", &record.generatedEnergy);
    tree->SetBranchAddress("GeneratedUx", &generatedUx);
    tree->SetBranchAddress("GeneratedUy", &generatedUy);
    tree->SetBranchAddress("GeneratedUz", &generatedUz);
    tree->SetBranchAddress("GeneratedTheta_rad", &record.generatedTheta);
    tree->SetBranchAddress("GeneratedPhi_rad", &generatedPhi);
    tree->SetBranchAddress("GeneratedSinTheta", &record.generatedSinTheta);
    tree->SetBranchAddress("GeneratedAbsUx", &generatedAbsUx);
    tree->SetBranchAddress("GasEntered", &record.gasEntered);
    tree->SetBranchAddress("EntryEnergy_keV", &record.entryEnergy);
    tree->SetBranchAddress("DepositedEnergy_keV", &record.depositedEnergy);
    tree->SetBranchAddress("Ux", &record.ux);
    tree->SetBranchAddress("Uy", &record.uy);
    tree->SetBranchAddress("Uz", &record.uz);
    tree->SetBranchAddress("Theta_rad", &record.theta);
    tree->SetBranchAddress("Phi_rad", &record.phi);
    tree->SetBranchAddress("SinTheta", &record.sinTheta);
    tree->SetBranchAddress("AbsUx", &absUx);
    tree->SetBranchAddress("FullyContained", &record.fullyContained);
    tree->SetBranchAddress("Selected", &record.selected);

    constexpr std::array<double, 6> edges = {{10., 15., 20., 30., 40., 60.}};
    constexpr std::array<Stage, 4> stages = {{
      Stage::Generated, Stage::GasEntry, Stage::Contained, Stage::Selected
    }};
    std::array<std::array<Accumulator, 5>, 4> accumulators;
    std::array<std::array<Summary, 5>, 4> summaries;
    long long angleConsistencyWarnings = 0;
    long long primaryBetaRows = 0;

    for (long long entry = 0; entry < tree->GetEntries(); ++entry) {
      tree->GetEntry(entry);
      // In this single-thread Sr-90 chain, primary ion TrackID=1 is Sr-90 and
      // TrackID=2 is Y-90.  Other RadioactiveDecay electrons are rare internal
      // conversion/de-excitation products and are not primary betas.
      if (record.parentID != 1 && record.parentID != 2) continue;
      ++primaryBetaRows;
      const auto generatedAngles =
        Angles({generatedUx, generatedUy, generatedUz});
      if (!NearlyEqual(generatedAngles.theta, record.generatedTheta, 1.e-10) ||
          !NearlyEqual(generatedAngles.phi, generatedPhi, 1.e-10) ||
          !NearlyEqual(generatedAngles.sinTheta, record.generatedSinTheta, 1.e-10) ||
          !NearlyEqual(generatedAngles.absUx, generatedAbsUx, 1.e-10)) {
        ++angleConsistencyWarnings;
      }
      if (record.gasEntered) {
        const auto entryAngles = Angles({record.ux, record.uy, record.uz});
        if (!NearlyEqual(entryAngles.theta, record.theta, 1.e-10) ||
            !NearlyEqual(entryAngles.phi, record.phi, 1.e-10) ||
            !NearlyEqual(entryAngles.sinTheta, record.sinTheta, 1.e-10) ||
            !NearlyEqual(entryAngles.absUx, absUx, 1.e-10)) {
          ++angleConsistencyWarnings;
        }
      }

      for (std::size_t stageIndex = 0; stageIndex < stages.size(); ++stageIndex) {
        const Stage stage = stages[stageIndex];
        if (!Passes(stage, record)) continue;
        const double energy = stage == Stage::Generated
          ? record.generatedEnergy : record.depositedEnergy;
        const double theta = stage == Stage::Generated
          ? record.generatedTheta : record.theta;
        const double sinTheta = stage == Stage::Generated
          ? record.generatedSinTheta : record.sinTheta;
        for (std::size_t bin = 0; bin < edges.size() - 1; ++bin) {
          if (!InOpenBin(energy, edges[bin], edges[bin + 1])) continue;
          auto& accumulator = accumulators[stageIndex][bin];
          accumulator.theta.push_back(theta);
          accumulator.sinTheta.push_back(sinTheta);
          accumulator.events.insert(record.eventID);
        }
      }
    }

    for (std::size_t stage = 0; stage < stages.size(); ++stage) {
      for (std::size_t bin = 0; bin < edges.size() - 1; ++bin) {
        const auto& accumulator = accumulators[stage][bin];
        summaries[stage][bin] =
          Summarize(accumulator.theta, accumulator.sinTheta, accumulator.events);
      }
    }

    std::ofstream csv(csvPath);
    if (!csv) throw std::runtime_error("cannot create CSV output");
    csv << "energy_low_keV,energy_high_keV,stage,energy_definition,"
           "electron_count,event_count,theta_median_deg,theta_p16_deg,"
           "theta_p84_deg,mean_sin_theta,rms_sin_theta,"
           "frac_abs_theta_minus_90_lt_10deg,"
           "frac_abs_theta_minus_90_lt_20deg,"
           "frac_abs_theta_minus_90_lt_30deg,"
           "frac_sin_theta_lt_0p25,frac_sin_theta_lt_0p50,"
           "frac_sin_theta_lt_0p75\n";
    csv << std::setprecision(12);
    for (std::size_t stage = 0; stage < stages.size(); ++stage) {
      for (std::size_t bin = 0; bin < edges.size() - 1; ++bin) {
        const auto& value = summaries[stage][bin];
        csv << edges[bin] << ',' << edges[bin + 1] << ','
            << StageName(stages[stage]) << ',' << EnergyDefinition(stages[stage])
            << ',' << value.electronCount << ',' << value.eventCount << ',';
        WriteNumber(csv, value.thetaMedianDeg); csv << ',';
        WriteNumber(csv, value.thetaP16Deg); csv << ',';
        WriteNumber(csv, value.thetaP84Deg); csv << ',';
        WriteNumber(csv, value.meanSinTheta); csv << ',';
        WriteNumber(csv, value.rmsSinTheta); csv << ',';
        WriteNumber(csv, value.fracPlane10); csv << ',';
        WriteNumber(csv, value.fracPlane20); csv << ',';
        WriteNumber(csv, value.fracPlane30); csv << ',';
        WriteNumber(csv, value.fracSin025); csv << ',';
        WriteNumber(csv, value.fracSin050); csv << ',';
        WriteNumber(csv, value.fracSin075); csv << '\n';
      }
    }
    csv.close();

    TFile output(rootPath.c_str(), "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("cannot create ROOT output");
    output.cd();
    auto* acceptance = tree->CopyTree("ParentID==1 || ParentID==2");
    acceptance->SetName("Acceptance");
    acceptance->SetTitle("One row per beta electron created by radioactive decay");
    acceptance->Write();

    TNamed definitions(
      "definitions",
      "analysis coordinates (x,y,z)=(G4 y,G4 x,G4 z); "
      "theta=acos(ux); phi=atan2(uy,uz); energy bins are strict open intervals; "
      "selected=FullyContained&&RecoAvailable; downstream bins use "
      "DepositedEnergy_keV as in the manuscript analysis; primary beta filter "
      "is ParentID==1 (Sr90) or ParentID==2 (Y90)");
    definitions.Write();

    double energyLow = 0.;
    double energyHigh = 0.;
    std::string stageName;
    std::string energyDefinition;
    long long electronCount = 0;
    long long eventCount = 0;
    Summary summaryRow;
    TTree summaryTree("Summary", "Acceptance summary by stage and energy bin");
    summaryTree.Branch("EnergyLow_keV", &energyLow);
    summaryTree.Branch("EnergyHigh_keV", &energyHigh);
    summaryTree.Branch("Stage", &stageName);
    summaryTree.Branch("EnergyDefinition", &energyDefinition);
    summaryTree.Branch("ElectronCount", &electronCount);
    summaryTree.Branch("EventCount", &eventCount);
    summaryTree.Branch("ThetaMedian_deg", &summaryRow.thetaMedianDeg);
    summaryTree.Branch("ThetaP16_deg", &summaryRow.thetaP16Deg);
    summaryTree.Branch("ThetaP84_deg", &summaryRow.thetaP84Deg);
    summaryTree.Branch("MeanSinTheta", &summaryRow.meanSinTheta);
    summaryTree.Branch("RmsSinTheta", &summaryRow.rmsSinTheta);
    summaryTree.Branch("FracPlane10", &summaryRow.fracPlane10);
    summaryTree.Branch("FracPlane20", &summaryRow.fracPlane20);
    summaryTree.Branch("FracPlane30", &summaryRow.fracPlane30);
    summaryTree.Branch("FracSin025", &summaryRow.fracSin025);
    summaryTree.Branch("FracSin050", &summaryRow.fracSin050);
    summaryTree.Branch("FracSin075", &summaryRow.fracSin075);
    for (std::size_t stage = 0; stage < stages.size(); ++stage) {
      for (std::size_t bin = 0; bin < edges.size() - 1; ++bin) {
        energyLow = edges[bin];
        energyHigh = edges[bin + 1];
        stageName = StageName(stages[stage]);
        energyDefinition = EnergyDefinition(stages[stage]);
        summaryRow = summaries[stage][bin];
        electronCount = summaryRow.electronCount;
        eventCount = summaryRow.eventCount;
        summaryTree.Fill();
      }
    }
    summaryTree.Write();

    auto makeCountHistogram = [&](const char* name, std::size_t stage) {
      auto* histogram = new TH1D(name, name, 5, 0., 5.);
      for (std::size_t bin = 0; bin < 5; ++bin) {
        histogram->SetBinContent(bin + 1, summaries[stage][bin].electronCount);
        histogram->GetXaxis()->SetBinLabel(
          bin + 1, (std::to_string(static_cast<int>(edges[bin])) + "-" +
                    std::to_string(static_cast<int>(edges[bin + 1]))).c_str());
      }
      histogram->GetXaxis()->SetTitle("Energy bin (keV)");
      histogram->GetYaxis()->SetTitle("Beta-electron records");
      return histogram;
    };

    auto* hGenerated = makeCountHistogram("generated_counts", 0);
    auto* hEntry = makeCountHistogram("gas_entry_counts", 1);
    auto* hContained = makeCountHistogram("contained_counts", 2);
    auto* hSelected = makeCountHistogram("selected_counts", 3);
    hGenerated->SetLineColor(kBlack);
    hEntry->SetLineColor(kBlue + 1);
    hContained->SetLineColor(kOrange + 7);
    hSelected->SetLineColor(kRed + 1);
    for (auto* histogram : {hGenerated, hEntry, hContained, hSelected}) {
      histogram->SetLineWidth(2);
      histogram->Write();
    }

    std::array<double, 5> centers;
    std::array<double, 5> halfWidths;
    std::array<double, 5> medians;
    std::array<double, 5> lowerErrors;
    std::array<double, 5> upperErrors;
    std::array<double, 5> meanSin;
    std::array<double, 5> rmsSin;
    std::array<double, 5> zeros = {{0.,0.,0.,0.,0.}};
    for (std::size_t bin = 0; bin < 5; ++bin) {
      centers[bin] = 0.5*(edges[bin] + edges[bin + 1]);
      halfWidths[bin] = 0.5*(edges[bin + 1] - edges[bin]);
      const auto& selected = summaries[3][bin];
      medians[bin] = selected.thetaMedianDeg;
      lowerErrors[bin] = selected.thetaMedianDeg - selected.thetaP16Deg;
      upperErrors[bin] = selected.thetaP84Deg - selected.thetaMedianDeg;
      meanSin[bin] = selected.meanSinTheta;
      rmsSin[bin] = selected.rmsSinTheta;
    }

    TGraphAsymmErrors thetaGraph(
      5, centers.data(), medians.data(), halfWidths.data(), halfWidths.data(),
      lowerErrors.data(), upperErrors.data());
    thetaGraph.SetName("selected_theta_quantiles");
    thetaGraph.SetTitle("Selected gas-entry polar angle;Deposited energy (keV);#theta (deg)");
    thetaGraph.SetMarkerStyle(20);
    thetaGraph.Write();

    TGraphErrors sinGraph(
      5, centers.data(), meanSin.data(), halfWidths.data(), rmsSin.data());
    sinGraph.SetName("selected_sin_theta");
    sinGraph.SetTitle("Selected gas-entry sin(#theta);Deposited energy (keV);sin(#theta)");
    sinGraph.SetMarkerStyle(20);
    sinGraph.Write();

    auto makeFractionGraph = [&](const char* name, const char* title,
                                 double Summary::*member, int color,
                                 int marker) {
      std::array<double, 5> values;
      for (std::size_t bin = 0; bin < 5; ++bin) {
        values[bin] = summaries[3][bin].*member;
      }
      auto* graph = new TGraphErrors(
        5, centers.data(), values.data(), halfWidths.data(), zeros.data());
      graph->SetName(name);
      graph->SetTitle(title);
      graph->SetLineColor(color);
      graph->SetMarkerColor(color);
      graph->SetMarkerStyle(marker);
      graph->SetLineWidth(2);
      graph->Write();
      return graph;
    };

    auto* plane10 = makeFractionGraph(
      "selected_frac_plane10", "|#theta-90#circ|<10#circ",
      &Summary::fracPlane10, kBlue + 1, 20);
    auto* plane20 = makeFractionGraph(
      "selected_frac_plane20", "|#theta-90#circ|<20#circ",
      &Summary::fracPlane20, kOrange + 7, 21);
    auto* plane30 = makeFractionGraph(
      "selected_frac_plane30", "|#theta-90#circ|<30#circ",
      &Summary::fracPlane30, kRed + 1, 22);
    auto* sin025 = makeFractionGraph(
      "selected_frac_sin025", "sin(#theta)<0.25",
      &Summary::fracSin025, kBlue + 1, 20);
    auto* sin050 = makeFractionGraph(
      "selected_frac_sin050", "sin(#theta)<0.50",
      &Summary::fracSin050, kOrange + 7, 21);
    auto* sin075 = makeFractionGraph(
      "selected_frac_sin075", "sin(#theta)<0.75",
      &Summary::fracSin075, kRed + 1, 22);

    gStyle->SetOptStat(0);
    TCanvas canvas("acceptance_canvas", "Sr90 polar-angle acceptance", 1800, 1100);
    canvas.Divide(3, 2);
    for (int pad = 1; pad <= 6; ++pad) {
      canvas.cd(pad);
      gPad->SetLeftMargin(0.14);
      gPad->SetRightMargin(0.04);
      gPad->SetBottomMargin(0.13);
      gPad->SetTopMargin(0.10);
    }
    canvas.cd(1);
    gPad->SetLogy();
    hGenerated->SetTitle("Acceptance counts;Energy bin (keV);Beta-electron records");
    hGenerated->Draw("HIST");
    hEntry->Draw("HIST SAME");
    hContained->Draw("HIST SAME");
    hSelected->Draw("HIST SAME");
    {
      TLegend legend(0.55, 0.64, 0.88, 0.88);
      legend.AddEntry(hGenerated, "generated", "l");
      legend.AddEntry(hEntry, "gas entry", "l");
      legend.AddEntry(hContained, "contained", "l");
      legend.AddEntry(hSelected, "selected", "l");
      legend.Draw();
      canvas.cd(2);
      thetaGraph.Draw("AP");
      canvas.cd(3);
      sinGraph.Draw("AP");
      canvas.cd(4);
      plane10->SetTitle("Near-plane fractions;Deposited energy (keV);Fraction");
      plane10->GetYaxis()->SetRangeUser(0., 1.05);
      plane10->Draw("APL");
      plane20->SetTitle("");
      plane30->SetTitle("");
      plane20->Draw("PL SAME");
      plane30->Draw("PL SAME");
      TLegend planeLegend(0.48, 0.18, 0.88, 0.40);
      planeLegend.AddEntry(plane10, "|#theta-90#circ|<10#circ", "lp");
      planeLegend.AddEntry(plane20, "|#theta-90#circ|<20#circ", "lp");
      planeLegend.AddEntry(plane30, "|#theta-90#circ|<30#circ", "lp");
      planeLegend.Draw();
      canvas.cd(5);
      sin025->SetTitle("Low-sin(#theta) fractions;Deposited energy (keV);Fraction");
      sin025->GetYaxis()->SetRangeUser(0., 1.05);
      sin025->Draw("APL");
      sin050->SetTitle("");
      sin075->SetTitle("");
      sin050->Draw("PL SAME");
      sin075->Draw("PL SAME");
      TLegend sinLegend(0.55, 0.18, 0.88, 0.40);
      sinLegend.AddEntry(sin025, "sin(#theta)<0.25", "lp");
      sinLegend.AddEntry(sin050, "sin(#theta)<0.50", "lp");
      sinLegend.AddEntry(sin075, "sin(#theta)<0.75", "lp");
      sinLegend.Draw();
      canvas.cd(6);
      TH1D selectedTheta(
        "selected_theta_all_bins",
        "Selected gas-entry #theta, 10-60 keV;#theta (deg);Beta-electron records",
        90, 0., 180.);
      for (std::size_t bin = 0; bin < 5; ++bin) {
        for (double theta : accumulators[3][bin].theta) {
          selectedTheta.Fill(theta*kRadToDeg);
        }
      }
      selectedTheta.Draw();
      selectedTheta.Write();
      canvas.cd();
      canvas.Write();
      canvas.SaveAs(pngPath.c_str());
    }

    output.Close();
    input.Close();

    std::cout << "INPUT_ACCEPTANCE_ROWS " << inputRows << '\n';
    std::cout << "PRIMARY_BETA_ROWS " << primaryBetaRows << '\n';
    std::cout << "ANGLE_CONSISTENCY_WARNINGS " << angleConsistencyWarnings << '\n';
    for (std::size_t bin = 0; bin < 5; ++bin) {
      std::cout << "BIN " << edges[bin] << " " << edges[bin + 1];
      for (std::size_t stage = 0; stage < stages.size(); ++stage) {
        std::cout << " " << StageName(stages[stage]) << "_electrons="
                  << summaries[stage][bin].electronCount
                  << " " << StageName(stages[stage]) << "_events="
                  << summaries[stage][bin].eventCount;
      }
      std::cout << '\n';
    }
    std::cout << "WROTE " << rootPath << '\n'
              << "WROTE " << csvPath << '\n'
              << "WROTE " << pngPath << '\n';
    return angleConsistencyWarnings == 0 ? 0 : 3;
  } catch (const std::exception& error) {
    std::cerr << "ERROR " << error.what() << '\n';
    return 1;
  }
}
