using System;

namespace MeanlineSolver.Core
{
    public static class VelocityTriangle
    {
        /// <summary>
        /// Compute absolute velocity magnitude C
        /// </summary>
        public static double ComputeAbsoluteVelocity(double Ca, double Ct)
        {
            return Math.Sqrt(Ca * Ca + Ct * Ct);
        }

        /// <summary>
        /// Compute relative velocity components W = C - U (vector)
        /// </summary>
        public static (double Wa, double Wt, double W) ComputeRelativeVelocity(double Ca, double Ct, double U)
        {
            double Wa = Ca;
            double Wt = Ct - U;
            double W = Math.Sqrt(Wa * Wa + Wt * Wt);

            return (Wa, Wt, W);
        }

        /// <summary>
        /// Compute flow angle α = atan(Ct/Ca)
        /// </summary>
        public static double ComputeAlpha(double Ca, double Ct)
        {
            return Math.Atan2(Ct, Ca);   // radians
        }

        /// <summary>
        /// Compute relative flow angle β = atan(Wt/Wa)
        /// </summary>
        public static double ComputeBeta(double Wa, double Wt)
        {
            return Math.Atan2(Wt, Wa);
        }

        /// <summary>
        /// Compute rotor-exit tangential velocity Ct2 from loading coefficient ψ.
        ///
        /// Simplified (constant mean radius):
        ///   ψ = Δh0/U^2 and Δh0 = U (Ct2 - Ct1)  => Ct2 = Ct1 + ψ U
        ///
        /// Meanline notebook form (allowing rm change):
        ///   Ct2 = ψ U2 + (rm1/rm2) Ct1
        /// </summary>
        public static double ComputeCt2(double Ct1, double psi, double U)
        {
            return Ct1 + psi * U;
        }

        /// <summary>
        /// Compute Ct2 with mean-radius change correction:
        /// Ct2 = ψ U2 + (rm1/rm2) Ct1
        /// </summary>
        public static double ComputeCt2(double Ct1, double psi, double U2, double rm1, double rm2)
        {
            if (rm2 == 0.0)
                throw new ArgumentException("rm2 cannot be zero.", nameof(rm2));

            return psi * U2 + (rm1 / rm2) * Ct1;
        }

        /// <summary>
        /// Compute axial velocity at next station (using Axial Velocity Ratio)
        /// Ca2 = Ca1 * AVR
        ///</summary>
        public static double ComputeNextAxialVelocity(double Ca1, double axialVelocityRatio)
        {
            return Ca1 * axialVelocityRatio;
        }

        /// <summary>
        /// Compute Euler turbine/compressor work:
        /// Δh0 = U (Ct2 - Ct1)
        /// </summary>
        public static double ComputeDeltaH0(double U, double Ct2, double Ct1)
        {
            return U * (Ct2 - Ct1);
        }

        /// <summary>
        /// General Euler work: Δh0 = U2*Ct2 - U1*Ct1
        /// </summary>
        public static double ComputeDeltaH0(double U1, double U2, double Ct1, double Ct2)
        {
            return U2 * Ct2 - U1 * Ct1;
        }

        /// <summary>
        /// Compute flow coefficient φ = Ca / U
        /// </summary>
        public static double ComputePhi(double Ca, double U)
        {
            return Ca / U;
        }

        /// <summary>
        /// Compute loading coefficient ψ = Δh0 / U^2
        /// </summary>
        public static double ComputePsi(double deltaH0, double U)
        {
            return deltaH0 / (U * U);
        }
    }
}
