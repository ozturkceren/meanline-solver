using System;

namespace MeanlineSolver.Core
{
    /// <summary>
    /// Empirical correlations for incidence, deviation and losses.
    /// NOTE: These are simplified placeholder models.
    /// You can later replace them with detailed correlations from the PDF.
    /// </summary>
    public static class EmpiricalModels
    {
        // -----------------------------
        // INCIDENCE ANGLE (i)
        // -----------------------------
        //
        // Very simple model:
        // i = K_inc * (phi - phi_opt)
        // where phi_opt depends weakly on reaction and solidity.
        //
        public static double ComputeIncidence(
            double flowCoefficient,   // φ
            double reaction,          // Λ
            double solidity)          // σ
        {
            // "Optimum" flow coefficient (very rough guess)
            double phiOpt = 0.5 + 0.1 * (reaction - 0.5); // around 0.5 for Λ=0.5
            double Kinc = 5.0 * Math.PI / 180.0;          // 5 deg per unit φ difference

            double i = Kinc * (flowCoefficient - phiOpt);

            return i; // radians
        }

        // -----------------------------
        // DEVIATION ANGLE (δ)
        // -----------------------------
        //
        // Simple Lieblein-like trend:
        // δ ∝ (t/c) / σ  (more thickness and less solidity → more deviation)
        //
        public static double ComputeDeviation(
            double thicknessChordRatio,  // t/c
            double solidity)             // σ
        {
            // Base deviation ~ 4 deg for typical σ and t/c
            double baseDeg = 4.0;

            double factor = (thicknessChordRatio / 0.06) * (1.2 / solidity);
            double deltaDeg = baseDeg * factor;

            return deltaDeg * Math.PI / 180.0; // radians
        }

        // -----------------------------
        // PROFILE LOSS COEFFICIENT (Y_p)
        // -----------------------------
        //
        // Very simplified profile loss:
        // Y_p ∝ (t/c) * (1/σ)
        //
        public static double ComputeProfileLoss(
            double thicknessChordRatio,  // t/c
            double solidity)             // σ
        {
            // base profile loss at design conditions
            double baseLoss = 0.03; // typical 3% total pressure loss

            double factor = (thicknessChordRatio / 0.06) * (1.2 / solidity);

            return baseLoss * factor;  // dimensionless loss coefficient
        }

        // -----------------------------
        // SECONDARY LOSS (Y_s)
        // -----------------------------
        //
        // Very rough model:
        // Y_s ∝ (1 - hubTipRatio)
        //
        public static double ComputeSecondaryLoss(double hubTipRatio)
        {
            double baseLoss = 0.01; // 1% at typical aspect
            double factor = (1.0 - hubTipRatio) / 0.2; // more span → more secondary loss

            return baseLoss * Math.Max(factor, 0.0);
        }

        // -----------------------------
        // TIP CLEARANCE LOSS (Y_c)
        // -----------------------------
        //
        // Y_c ∝ (clearance / bladeHeight)
        //
        public static double ComputeClearanceLoss(
            double tipClearance,  // absolute clearance
            double bladeHeight)   // span at tip
        {
            if (bladeHeight <= 0.0)
                return 0.0;

            double clearanceRatio = tipClearance / bladeHeight;

            double baseLoss = 0.015; // 1.5% for typical clearance
            double factor = clearanceRatio / 0.01; // 1% span as reference

            return baseLoss * factor;
        }

        // -----------------------------
        // TOTAL LOSS COEFFICIENT (Y_total)
        // -----------------------------
        //
        // Simple sum of components:
        //
        public static double ComputeTotalLoss(
            double thicknessChordRatio,
            double solidity,
            double hubTipRatio,
            double tipClearance,
            double bladeHeight)
        {
            double Yp = ComputeProfileLoss(thicknessChordRatio, solidity);
            double Ys = ComputeSecondaryLoss(hubTipRatio);
            double Yc = ComputeClearanceLoss(tipClearance, bladeHeight);

            return Yp + Ys + Yc;
        }

        // -----------------------------
        // SIMPLE EFFICIENCY FROM LOSS
        // -----------------------------
        //
        // Very rough link between loss coefficient and efficiency.
        //
        public static double EstimateEfficiencyFromLoss(double totalLoss)
        {
            // Assume loss coefficient around 0.05 → efficiency ~ 0.9
            double eff = 1.0 - 0.5 * totalLoss;

            // clamp
            if (eff > 0.99) eff = 0.99;
            if (eff < 0.70) eff = 0.70;

            return eff;
        }
    }
}
