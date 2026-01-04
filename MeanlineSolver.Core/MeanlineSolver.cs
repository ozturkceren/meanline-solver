using System;
using System.Collections.Generic;

namespace MeanlineSolver.Core
{
    public class MeanlineSolverMain
    {
        private readonly StageSolver _stageSolver;

        public MeanlineSolverMain()
        {
            _stageSolver = new StageSolver();
        }

        public List<StageResult> Solve(MeanlineInputs inputs)
        {
            var results = new List<StageResult>();

            // Global inlet total conditions
            double T0_in = inputs.InletTotalTemperature;
            double P0_in = inputs.InletTotalPressure;

            // Starting radius ratios remain from inputs
            double hubTip = inputs.InletHubTipRatio;

            for (int i = 0; i < inputs.Stages.Count; i++)
            {
                var stageInput = inputs.Stages[i];

                // Update global inlet for this stage
                inputs.InletTotalTemperature = T0_in;
                inputs.InletTotalPressure = P0_in;
                inputs.InletHubTipRatio = hubTip;

                // Solve single stage
                var result = _stageSolver.SolveSingleStage(stageInput, inputs);
                results.Add(result);

                // --- Compute new inlet for next stage ---
                // Use station 3 values
                var st3 = result.Station3;

                // Total values propagate
                T0_in = st3.T0;
                P0_in = st3.P0;

                // Ca₁ next stage = Ca₃ previous stage
                // This is implicitly handled inside ComputeStation1 using φ, U

                // radius update (very simple meanline assumption)
                hubTip = st3.Height / (st3.Radius + st3.Height);

                // prepare for next iteration
            }

            return results;
        }
    }
}
