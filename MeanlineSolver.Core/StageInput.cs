namespace MeanlineSolver.Core
{
    /// <summary>
    /// Stage-wise inputs created strictly from the source INPUT list.
    /// Outer-loop variables:
    ///  - LoadingCoefficient (psi)   : adjusted by overall pressure ratio check
    ///  - StatorExitFlowAngleDeg     : adjusted by stage reaction check
    /// </summary>
    public sealed class StageInput
    {
        public int Index { get; init; }

        // From INPUT list
        public double Reaction { get; init; }            // target stage reaction R_target
        public double LoadingCoefficient { get; init; }  // stage loading ψ (outer-loop variable)

        public double FlowCoefficient { get; init; }     // inlet flow coefficient φ
        public double HubTipRatio { get; init; }         // inlet hub-tip ratio rh/rt
        public double AxialVelocityRatio { get; init; }  // AVR (single input)
        public double InletFlowAngleDeg { get; init; }   // α1 [deg]

        // Outer-loop variable for reaction check
        public double StatorExitFlowAngleDeg { get; init; } // α3 guess [deg]
    }
}
