using MeanlineSolver.Core;

public class StageResult
{
    public StationState Station1 { get; set; }
    public StationState Station2 { get; set; }
    public StationState Station3 { get; set; }

    public double DeltaH0 { get; set; }
    public double PressureRatio { get; set; }
    public double Efficiency { get; set; }

    public double FlowCoefficient { get; set; }
    public double LoadingCoefficient { get; set; }
    public double Reaction { get; set; }

    public double RotorInletAngle { get; set; }
    public double RotorExitAngle { get; set; }
    public double StatorInletAngle { get; set; }
    public double StatorExitAngle { get; set; }

    public double Incidence { get; set; }
    public double Deviation { get; set; }
    public double Solidity { get; set; }

    // ADD THESE (pdf output compatibility)
    public double ChordLengthRotor { get; set; }
    public double ChordLengthStator { get; set; }
    public double CamberAngleRotor { get; set; }
    public double CamberAngleStator { get; set; }
}
