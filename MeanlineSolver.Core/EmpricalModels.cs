using System;

namespace MeanlineSolver.Core
{
    /// <summary>
    /// Correlation layer constrained to the INPUT list shown in the source screenshot.
    /// The exact report correlations are not visible in the image,
    /// so mappings below are deterministic and only use provided inputs.
    /// </summary>
    public static class EmpiricalModels
    {
        // DF -> solidity mapping (deterministic)
        public static double ComputeSolidityFromDiffusionFactor(double diffusionFactor)
        {
            double df = Clamp(diffusionFactor, 0.2, 0.7);
            double sigma = 0.8 + (df - 0.2) * (2.0 - 0.8) / (0.7 - 0.2);
            return Clamp(sigma, 0.6, 2.2);
        }

        // i* (rad)
        // 3-arg overload (IterateAngles.cs bunu çağırıyor)
        public static double ComputeIncidence(double phi, double reactionTarget, double solidity)
        {
            // Sadece INPUT listesindeki büyüklükleri kullanır:
            // phi (Inlet Flow Coefficient), reactionTarget (Stage Reactions), solidity (Solidity output/iter)
            // Deterministik, bounded bir korelasyon (defter mantığıyla "i* hesapla" adımını sağlar)

            double ph = Clamp(phi, 0.2, 0.9);
            double r = Clamp(reactionTarget, 0.2, 0.8);
            double s = Clamp(solidity, 0.5, 2.5);

            // Basit ama stabil bir model:
            // phi ↑ => incidence biraz artar
            // reaction ↑ => incidence biraz azalır (tipik tasarım trendi)
            // solidity ↑ => incidence biraz azalır (daha yüksek σ genelde daha toleranslı)
            double deg =
                6.0 * (ph - 0.5)
              - 5.0 * (r - 0.5)
              - 2.0 * (s - 1.0);

            deg = Clamp(deg, -6.0, 6.0);
            return deg * Math.PI / 180.0;
        }

        // δ* (rad)
        public static double ComputeDeviation(double diffusionFactor, double thicknessChordRatio)
        {
            double df = Clamp(diffusionFactor, 0.2, 0.7);
            double tc = Clamp(thicknessChordRatio, 0.02, 0.12);

            double deg = 2.0 + 12.0 * (df - 0.3) + 25.0 * (tc - 0.06);
            deg = Clamp(deg, 0.0, 12.0);
            return deg * Math.PI / 180.0;
        }

        // Rotor loss Y (bounded)
        public static double ComputeLossY_Rotor(double diffusionFactor, double thicknessChordRatio, double tipClearance, double bladeHeight)
        {
            double df = Clamp(diffusionFactor, 0.0, 0.95);
            double tc = Clamp(thicknessChordRatio, 0.0, 0.2);
            double clr = (bladeHeight > 1e-9) ? Clamp(tipClearance / bladeHeight, 0.0, 0.2) : 0.0;

            double y = 0.06 + 0.18 * df + 0.25 * tc + 0.30 * clr;
            return Clamp(y, 0.0, 0.35);
        }

        // Stator loss Y (bounded)
        public static double ComputeLossY_Stator(double diffusionFactor, double thicknessChordRatio)
        {
            double df = Clamp(diffusionFactor, 0.0, 0.95);
            double tc = Clamp(thicknessChordRatio, 0.0, 0.2);
            double y = 0.05 + 0.20 * df + 0.22 * tc;
            return Clamp(y, 0.0, 0.30);
        }

        // Δs = -R ln(1 - Y)
        public static double ComputeEntropyRiseFromY(double lossY)
        {
            double y = Clamp(lossY, 0.0, 0.95);
            return -FlowProperties.R * Math.Log(Math.Max(1.0 - y, 1e-12));
        }

        private static double Clamp(double x, double lo, double hi)
        {
            if (x < lo) return lo;
            if (x > hi) return hi;
            return x;
        }
    }
}
