using System;
using System.Collections.Generic;

namespace MeanlineSolver.Core
{
    public sealed class MeanlineSolverMain
    {
        private readonly StageSolver _stageSolver = new StageSolver();

        public List<StageResult> Solve(MeanlineInputs inputs)
        {
            ValidateInputs(inputs);

            var psiBase = BuildPsiFromDistribution(inputs);

            double psiScale = 1.0;
            List<StageResult> best = new();

            for (int outer = 0; outer < 12; outer++)
            {
                var results = new List<StageResult>(inputs.NumberOfStages);

                double T0 = inputs.InletTotalTemperature;
                double P0 = inputs.InletTotalPressure;

                for (int i = 0; i < inputs.NumberOfStages; i++)
                {
                    var stg = new StageInput
                    {
                        Index = i + 1,
                        Reaction = inputs.StageReactions[i],
                        LoadingCoefficient = psiBase[i] * psiScale,
                        FlowCoefficient = inputs.InletFlowCoefficient,
                        HubTipRatio = inputs.InletHubTipRatio,
                        AxialVelocityRatio = inputs.AxialVelocityRatio,
                        InletFlowAngleDeg = inputs.InletFlowAngle
                    };

                    var res = _stageSolver.SolveSingleStage(stg, inputs, T0, P0);
                    results.Add(res);

                    T0 = res.Station3.T0;
                    P0 = res.Station3.P0;
                }

                double PRcalc = results[^1].Station3.P0 / results[0].Station1.P0;
                best = results;

                double err = (PRcalc - inputs.PressureRatio) / inputs.PressureRatio;
                if (Math.Abs(err) < 1e-4)
                    break;

                psiScale *= Math.Pow(inputs.PressureRatio / Math.Max(PRcalc, 1e-12), 0.45);
                psiScale = Math.Clamp(psiScale, 0.2, 5.0);
            }

            return best;
        }

        private static void ValidateInputs(MeanlineInputs i)
        {
            if (i.NumberOfStages <= 0)
                throw new ArgumentException("NumberOfStages must be > 0");
            if (i.StageReactions == null || i.StageReactions.Length != i.NumberOfStages)
                throw new ArgumentException("StageReactions must have length = NumberOfStages");
            if (i.StageLoadDistribution == null || i.StageLoadDistribution.Length != i.NumberOfStages)
                throw new ArgumentException("StageLoadDistribution must have length = NumberOfStages");
        }

        private static double[] BuildPsiFromDistribution(MeanlineInputs i)
        {
            double sum = 0.0;
            for (int k = 0; k < i.StageLoadDistribution.Length; k++) sum += i.StageLoadDistribution[k];

            var psi = new double[i.NumberOfStages];

            if (Math.Abs(sum - 1.0) < 1e-6)
            {
                double psiMean = 0.35; // deterministic baseline
                for (int k = 0; k < psi.Length; k++)
                    psi[k] = psiMean * i.NumberOfStages * i.StageLoadDistribution[k];
            }
            else
            {
                for (int k = 0; k < psi.Length; k++)
                    psi[k] = i.StageLoadDistribution[k];
            }

            return psi;
        }
    }
}
