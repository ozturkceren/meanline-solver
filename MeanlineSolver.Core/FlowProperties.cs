using System;

namespace MeanlineSolver.Core
{
    public static class FlowProperties
    {
        public const double Gamma = 1.4;
        public const double R = 287.0;
        public const double Cp = Gamma * R / (Gamma - 1.0);

        // s = Cp ln(T) - R ln(P)  (constant cancels in differences)
        public static double ComputeEntropy(double T, double P)
        {
            return Cp * Math.Log(Math.Max(T, 1e-12)) - R * Math.Log(Math.Max(P, 1e-12));
        }

        // P = exp((Cp ln(T) - s)/R)
        public static double ComputePressureFromEntropy(double T, double s)
        {
            return Math.Exp((Cp * Math.Log(Math.Max(T, 1e-12)) - s) / R);
        }

        public static double ComputeTotalTemperature(double T, double Mach)
        {
            return T * (1.0 + (Gamma - 1.0) / 2.0 * Mach * Mach);
        }

        public static double ComputeTotalPressure(double P, double T, double T0)
        {
            double ratio = T0 / T;
            return P * Math.Pow(ratio, Gamma / (Gamma - 1.0));
        }

        public static double ComputeStaticTemperature(double T0, double velocity)
        {
            return T0 - velocity * velocity / (2.0 * Cp);
        }

        public static double ComputeStaticPressure(double P0, double T, double T0)
        {
            double ratio = T / T0;
            return P0 * Math.Pow(ratio, Gamma / (Gamma - 1.0));
        }

        public static double ComputeDensity(double P, double T)
        {
            return P / (R * T);
        }

        public static double ComputeSpeedOfSound(double T)
        {
            return Math.Sqrt(Gamma * R * T);
        }

        public static double ComputeMach(double velocity, double temperature)
        {
            double a = ComputeSpeedOfSound(temperature);
            return velocity / a;
        }

        public static double ComputeEnthalpy(double T)
        {
            return Cp * T;
        }

        public static double ComputeTotalEnthalpy(double h, double velocity)
        {
            return h + 0.5 * velocity * velocity;
        }

        public static double ComputeIsentropicPressureRatio(double deltaH0, double T01)
        {
            double term = 1.0 + deltaH0 / (Cp * T01);
            return Math.Pow(term, Gamma / (Gamma - 1.0));
        }
    }
}
