using System;

namespace MeanlineSolver.Core
{
    public static class FlowProperties
    {
        // GAS CONSTANTS (You can move these to MeanlineInputs if needed)
        public const double Gamma = 1.4;      // heat capacity ratio
        public const double R = 287.0;        // gas constant for air (J/kg K)
        public const double Cp = Gamma * R / (Gamma - 1.0);

        // ------------------------------
        //   ISENTROPIC RELATIONS
        // ------------------------------

        /// <summary>
        /// T0 = T * (1 + (gamma-1)/2 * M^2)
        /// </summary>
        public static double ComputeTotalTemperature(double T, double Mach)
        {
            return T * (1.0 + (Gamma - 1.0) / 2.0 * Mach * Mach);
        }

        /// <summary>
        /// P0 = P * (T0/T)^(gamma/(gamma-1))
        /// </summary>
        public static double ComputeTotalPressure(double P, double T, double T0)
        {
            double ratio = T0 / T;
            return P * Math.Pow(ratio, Gamma / (Gamma - 1.0));
        }

        /// <summary>
        /// Static Temperature from total temperature and velocity:
        /// T = T0 - C^2 / (2*Cp)
        /// </summary>
        public static double ComputeStaticTemperature(double T0, double velocity)
        {
            return T0 - velocity * velocity / (2.0 * Cp);
        }

        /// <summary>
        /// Static pressure using isentropic relation:
        /// P = P0 * (T/T0)^(gamma/(gamma-1))
        /// </summary>
        public static double ComputeStaticPressure(double P0, double T, double T0)
        {
            double ratio = T / T0;
            return P0 * Math.Pow(ratio, Gamma / (Gamma - 1.0));
        }

        // ------------------------------
        //   VELOCITY, MACH, DENSITY
        // ------------------------------

        /// <summary>
        /// Density using ideal gas law: ρ = P / (R * T)
        /// </summary>
        public static double ComputeDensity(double P, double T)
        {
            return P / (R * T);
        }

        /// <summary>
        /// Speed of sound: a = sqrt(gamma * R * T)
        /// </summary>
        public static double ComputeSpeedOfSound(double T)
        {
            return Math.Sqrt(Gamma * R * T);
        }

        /// <summary>
        /// Mach number: M = C / a
        /// </summary>
        public static double ComputeMach(double velocity, double temperature)
        {
            double a = ComputeSpeedOfSound(temperature);
            return velocity / a;
        }

        // ------------------------------
        //   ENTHALPY & ENERGY
        // ------------------------------

        /// <summary>
        /// Static enthalpy: h = Cp * T
        /// </summary>
        public static double ComputeEnthalpy(double T)
        {
            return Cp * T;
        }

        /// <summary>
        /// Total enthalpy: h0 = h + C^2/2
        /// </summary>
        public static double ComputeTotalEnthalpy(double h, double velocity)
        {
            return h + 0.5 * velocity * velocity;
        }

        // ------------------------------
        //   ISENTROPIC PRESSURE RATIO
        // ------------------------------

        /// <summary>
        /// Isentropic total pressure ratio from Δh0:
        /// PR = (1 + Δh0/(Cp*T01))^(gamma/(gamma-1))
        /// </summary>
        public static double ComputeIsentropicPressureRatio(double deltaH0, double T01)
        {
            double term = 1.0 + deltaH0 / (Cp * T01);
            return Math.Pow(term, Gamma / (Gamma - 1.0));
        }
    }
}
