using System;

namespace MeanlineSolver.Core
{
    /// <summary>
    /// Iterative blade-metal-angle solver (very compact meanline model).
    ///
    /// Goals:
    /// - Convert flow angles (α/β) to blade metal angles by applying incidence i* and deviation δ*
    /// - Optionally iterate solidity (σ) if a design rule updates σ based on δ*
    ///
    /// Notes:
    /// - All angles are in radians inside this class.
    /// - Incidence / deviation correlations are in EmpiricalModels (placeholders).
    /// </summary>
    public static class IterateAngles
    {
        public sealed class BladeAngleResult
        {
            public double RotorInletMetal { get; init; }   // β1m
            public double RotorExitMetal { get; init; }    // β2m
            public double StatorInletMetal { get; init; }  // α2m
            public double StatorExitMetal { get; init; }   // α3m

            public double Incidence { get; init; }         // i* (rad)
            public double Deviation { get; init; }         // δ* (rad)
            public double Solidity { get; init; }          // σ (final)
        }

        /// <summary>
        /// Solve blade metal angles with a small fixed-point iteration on solidity.
        /// </summary>
        public static BladeAngleResult Solve(
            double phi,
            double reactionTarget,
            double thicknessChordRatio,
            double solidityInitial,
            double alpha1,
            double beta1,
            double alpha2,
            double beta2,
            double alpha3,
            int maxIter = 30,
            double sigmaTol = 1e-6)
        {
            double sigma = Math.Max(solidityInitial, 0.2);
            double incidence = 0.0;
            double deviation = 0.0;

            for (int k = 0; k < maxIter; k++)
            {
                double sigmaOld = sigma;

                incidence = EmpiricalModels.ComputeIncidence(phi, reactionTarget, sigma);
                deviation = EmpiricalModels.ComputeDeviation(thicknessChordRatio, sigma);

                // Metal angles (simple textbook convention)
                // Rotor metal angles are based on relative flow (β)
                double beta1m = beta1 - incidence;
                double beta2m = beta2 + deviation;

                // Stator metal angles are based on absolute flow (α)
                double alpha2m = alpha2;             // stator inlet metal ~ flow inlet (can be extended)
                double alpha3m = alpha3 + deviation; // simple deviation application

                // Simple solidity update rule (placeholder): higher deviation -> slightly higher σ
                sigma = sigmaOld * (1.0 + 0.5 * Math.Abs(deviation));
                if (Math.Abs(sigma - sigmaOld) < sigmaTol)
                {
                    return new BladeAngleResult
                    {
                        RotorInletMetal = beta1m,
                        RotorExitMetal = beta2m,
                        StatorInletMetal = alpha2m,
                        StatorExitMetal = alpha3m,
                        Incidence = incidence,
                        Deviation = deviation,
                        Solidity = sigma
                    };
                }
            }

            // If not converged, return last iteration
            return new BladeAngleResult
            {
                RotorInletMetal = beta1 - incidence,
                RotorExitMetal = beta2 + deviation,
                StatorInletMetal = alpha2,
                StatorExitMetal = alpha3 + deviation,
                Incidence = incidence,
                Deviation = deviation,
                Solidity = sigma
            };
        }
    }
}
