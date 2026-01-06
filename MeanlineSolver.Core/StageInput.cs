namespace MeanlineSolver.Core
{
    /// <summary>
    /// Stage-wise meanline inputs coming strictly from the source INPUT list.
    ///
    /// IMPORTANT:
    /// - No extra stage-only unknowns are introduced.
    /// - Reaction and Loading are taken from MeanlineInputs arrays.
    /// - FlowCoefficient is inlet flow coefficient (same definition as source input list).
    /// - HubTipRatio is taken as inlet hub-tip ratio and kept constant across stations.
    /// - AxialVelocityRatio is taken from the single source input "Axial Velocity Ratio".
    /// </summary>
    public sealed class StageInput
    {
        public int Index { get; init; }

        // Source inputs
        public double Reaction { get; init; }           // stage reaction R
        public double LoadingCoefficient { get; init; } // stage loading ψ

        public double FlowCoefficient { get; init; }    // inlet flow coefficient φ
        public double HubTipRatio { get; init; }        // inlet hub-tip ratio (rh/rt)
        public double AxialVelocityRatio { get; init; } // AVR (applied rotor & stator)
        public double InletFlowAngleDeg { get; init; }  // α1 [deg]
    }
}
