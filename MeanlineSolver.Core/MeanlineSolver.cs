using System;
using System.Collections.Generic;

namespace MeanlineSolver.Core
{
    public sealed class MeanlineSolverMain
    {
        private readonly StageSolver _stageSolver = new StageSolver();

        // Outer-loop controls (deterministic constants, not user inputs)
        private const int MaxPsiIterations = 20;
        private const int MaxAlphaIterations = 30;

        private const double ReactionTol = 1e-4;
        private const double PressureRatioTol = 1e-4;

        private const double AlphaGainDeg = 8.0;   // α3 update gain (deg per reaction error)
        private const double PsiGain = 0.45;       // ψ scaling exponent

        public List<StageResult> Solve(MeanlineInputs inputs)
        {
            ValidateInputs(inputs);

            // Base psi distribution from input (StageLoadDistribution)
            double[] psiBase = BuildPsiFromDistribution(inputs);

            // α3 initial guesses: use inlet flow angle (available input), deterministic, no new unknown input
            double[] alpha3Deg = new double[inputs.NumberOfStages];
            for (int i = 0; i < alpha3Deg.Length; i++)
                alpha3Deg[i] = inputs.InletFlowAngle; // initial guess

            double psiScale = 1.0;

            List<StageResult> best = new();

            for (int psiIter = 0; psiIter < MaxPsiIterations; psiIter++)
            {
                // ==========================================================
                // INNER LOOP: Adjust Stator Exit Flow Angle(s) alpha3 to match Reaction(s)
                // ==========================================================
                List<StageResult> resultsAfterAlpha = new();

                for (int aIter = 0; aIter < MaxAlphaIterations; aIter++)
                {
                    var stageResults = new List<StageResult>(inputs.NumberOfStages);

                    double T0 = inputs.InletTotalTemperature;
                    double P0 = inputs.InletTotalPressure;

                    double maxReactionAbsErr = 0.0;

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
                            InletFlowAngleDeg = inputs.InletFlowAngle,

                            StatorExitFlowAngleDeg = alpha3Deg[i]
                        };

                        double Rcalc;
                        StageResult res = _stageSolver.SolveSingleStage(stg, inputs, T0, P0, out Rcalc);
                        stageResults.Add(res);

                        // propagate totals to next stage
                        T0 = res.Station3.T0;
                        P0 = res.Station3.P0;

                        // reaction error for this stage
                        double err = stg.Reaction - Rcalc;
                        maxReactionAbsErr = Math.Max(maxReactionAbsErr, Math.Abs(err));

                        // Update alpha3 for this stage (ONLY alpha3 changes here)
                        // alpha3_new = alpha3_old + K * (Rtarget - Rcalc)
                        alpha3Deg[i] = Clamp(alpha3Deg[i] + AlphaGainDeg * err, -75.0, 75.0);
                    }

                    resultsAfterAlpha = stageResults;

                    if (maxReactionAbsErr < ReactionTol)
                        break;
                }

                best = resultsAfterAlpha;

                // ==========================================================
                // OUTER CHECK: Overall Pressure Ratio Check (update ONLY psiScale)
                // ==========================================================
                double PRcalc = best[^1].Station3.P0 / best[0].Station1.P0;
                double PRtarget = inputs.PressureRatio;

                double prErr = (PRcalc - PRtarget) / Math.Max(PRtarget, 1e-12);

                if (Math.Abs(prErr) < PressureRatioTol)
                    break;

                // Update psiScale (ONLY psi changes here)
                psiScale *= Math.Pow(PRtarget / Math.Max(PRcalc, 1e-12), PsiGain);
                psiScale = Clamp(psiScale, 0.2, 5.0);
            }

            return best;
        }

        private static void ValidateInputs(MeanlineInputs i)
        {
            if (i.NumberOfStages <= 0)
                throw new ArgumentException("NumberOfStages must be > 0");

            if (i.StageReactions == null || i.StageReactions.Length != i.NumberOfStages)
                throw new ArgumentException("StageReactions length must equal NumberOfStages");

            if (i.StageLoadDistribution == null || i.StageLoadDistribution.Length != i.NumberOfStages)
                throw new ArgumentException("StageLoadDistribution length must equal NumberOfStages");
        }

        private static double[] BuildPsiFromDistribution(MeanlineInputs i)
        {
            double sum = 0.0;
            for (int k = 0; k < i.StageLoadDistribution.Length; k++)
                sum += i.StageLoadDistribution[k];

            var psi = new double[i.NumberOfStages];

            // If distribution sums to 1 -> use it as weights around a deterministic baseline
            if (Math.Abs(sum - 1.0) < 1e-9)
            {
                // baseline stage loading (deterministic numerical seed)
                // NOTE: this is NOT a new user input; it is only the initial numeric guess.
                double psiMean = 0.35;

                for (int k = 0; k < psi.Length; k++)
                    psi[k] = psiMean * i.NumberOfStages * i.StageLoadDistribution[k];
            }
            else
            {
                // Otherwise treat values as absolute psi guesses
                for (int k = 0; k < psi.Length; k++)
                    psi[k] = i.StageLoadDistribution[k];
            }

            return psi;
        }

        private static double Clamp(double x, double lo, double hi)
        {
            if (x < lo) return lo;
            if (x > hi) return hi;
            return x;
        }
    }
}
