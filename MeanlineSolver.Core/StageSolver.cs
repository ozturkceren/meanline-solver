using System;

namespace MeanlineSolver.Core
{
    public sealed class StageSolver
    {
        private const int InletMaxIter = 60;
        private const int GeoMaxIter = 60;

        public StageResult SolveSingleStage(StageInput stg, MeanlineInputs g, double T01, double P01)
        {
            var s1 = ComputeStation1(stg, g, T01, P01);
            var s2 = ComputeStation2(stg, g, s1);
            var s3 = ComputeStation3(stg, g, s2);

            double PR_stage = s3.P0 / s1.P0;
            double deltaH0 = s2.H0 - s1.H0;

            var blade = ComputeBladeAngles(stg, g, s1, s2, s3);

            return new StageResult
            {
                Station1 = s1,
                Station2 = s2,
                Station3 = s3,

                DeltaH0 = deltaH0,
                PressureRatio = PR_stage,
                Efficiency = 1.0, // PR already computed from losses/entropy

                FlowCoefficient = stg.FlowCoefficient,
                LoadingCoefficient = stg.LoadingCoefficient,
                Reaction = stg.Reaction,

                RotorInletAngle = blade.RotorInletMetal,
                RotorExitAngle = blade.RotorExitMetal,
                StatorInletAngle = blade.StatorInletMetal,
                StatorExitAngle = blade.StatorExitMetal,

                Incidence = blade.Incidence,
                Deviation = blade.Deviation,
                Solidity = blade.Solidity,

                ChordLengthRotor = blade.ChordRotor,
                ChordLengthStator = blade.ChordStator,
                CamberAngleRotor = blade.CamberRotor,
                CamberAngleStator = blade.CamberStator,
            };
        }

        // =========================
        // STATION 1 : Inlet Iteration
        // =========================
        private StationState ComputeStation1(StageInput stg, MeanlineInputs g, double T01, double P01)
        {
            double omega = g.RotationalSpeed * 2.0 * Math.PI / 60.0;
            double ht = stg.HubTipRatio;
            double alpha1 = stg.InletFlowAngleDeg * Math.PI / 180.0;

            double rm = 0.25; // just an initial guess

            StationState s = new StationState { T0 = T01, P0 = P01 };

            for (int iter = 0; iter < InletMaxIter; iter++)
            {
                double U = omega * rm;
                double Ca = stg.FlowCoefficient * U;   // Ca = φU
                double Ct = Ca * Math.Tan(alpha1);     // Ct = Ca tan α1
                double C = Math.Sqrt(Ca * Ca + Ct * Ct);

                double T = FlowProperties.ComputeStaticTemperature(T01, C);
                double P = FlowProperties.ComputeStaticPressure(P01, T, T01);
                double rho = FlowProperties.ComputeDensity(P, T);

                double A = g.MassFlowRate / (rho * Ca * g.BlockageFactor);

                // A = 4π rm^2 (1-ht)/(1+ht)
                double rm_new = Math.Sqrt(A * (1.0 + ht) / (4.0 * Math.PI * (1.0 - ht)));

                if (Math.Abs(rm_new - rm) / Math.Max(rm, 1e-9) < 1e-9)
                {
                    rm = rm_new;
                    FillStationGeometry(ref s, rm, ht, A);
                    FillKinematicsAndThermo(ref s, omega, Ca, Ct, alpha1, T01, P01);
                    return s;
                }

                rm = 0.5 * (rm + rm_new);
            }

            // fallback
            double U_last = omega * rm;
            double Ca_last = stg.FlowCoefficient * U_last;
            double Ct_last = Ca_last * Math.Tan(alpha1);
            double C_last = Math.Sqrt(Ca_last * Ca_last + Ct_last * Ct_last);
            double T_last = FlowProperties.ComputeStaticTemperature(T01, C_last);
            double P_last = FlowProperties.ComputeStaticPressure(P01, T_last, T01);
            double rho_last = FlowProperties.ComputeDensity(P_last, T_last);
            double A_last = g.MassFlowRate / (rho_last * Ca_last * g.BlockageFactor);

            FillStationGeometry(ref s, rm, ht, A_last);
            FillKinematicsAndThermo(ref s, omega, Ca_last, Ct_last, alpha1, T01, P01);
            return s;
        }

        // =========================
        // STATION 2 : Rotor Exit
        // =========================
        private StationState ComputeStation2(StageInput stg, MeanlineInputs g, StationState s1)
        {
            double omega = g.RotationalSpeed * 2.0 * Math.PI / 60.0;
            double ht = stg.HubTipRatio;

            double rm2 = s1.Radius;
            double s2_static = s1.S;

            for (int iter = 0; iter < GeoMaxIter; iter++)
            {
                double U1 = s1.U;
                double U2 = omega * rm2;

                double Ca2 = stg.AxialVelocityRatio * s1.Ca;

                // Ct2 = ψ U2 + (rm1/rm2) Ct1
                double Ct2 = stg.LoadingCoefficient * U2 + (s1.Radius / rm2) * s1.Ct;

                double C2 = Math.Sqrt(Ca2 * Ca2 + Ct2 * Ct2);
                double alpha2 = Math.Atan2(Ct2, Ca2);

                double Wt2 = U2 - Ct2;
                double W2 = Math.Sqrt(Ca2 * Ca2 + Wt2 * Wt2);
                double beta2 = Math.Atan2(Wt2, Ca2);

                // Rothalpy I = h01 - U1 Ct1
                double I = s1.H0 - U1 * s1.Ct;

                // h2 = I - W2^2/2 + U2^2/2
                double h2 = I - 0.5 * W2 * W2 + 0.5 * U2 * U2;
                double T2 = h2 / FlowProperties.Cp;

                double h02 = h2 + 0.5 * C2 * C2;
                double T02 = h02 / FlowProperties.Cp;

                // provisional height for clearance loss normalization
                double P2_iso = s1.P0 * Math.Pow(T2 / s1.T0, FlowProperties.Gamma / (FlowProperties.Gamma - 1.0));
                double rho2_iso = FlowProperties.ComputeDensity(P2_iso, T2);
                double A2_iso = g.MassFlowRate / (rho2_iso * Ca2 * g.BlockageFactor);
                double h2_geo_iso = A2_iso / (2.0 * Math.PI * rm2);

                double Yrot = EmpiricalModels.ComputeLossY_Rotor(g.DiffusionFactor, g.ThicknessChordRatio, g.TipClearance, h2_geo_iso);
                double dS12 = EmpiricalModels.ComputeEntropyRiseFromY(Yrot);
                s2_static = s1.S + dS12;

                double P2 = FlowProperties.ComputePressureFromEntropy(T2, s2_static);
                double rho2 = FlowProperties.ComputeDensity(P2, T2);

                double A2 = g.MassFlowRate / (rho2 * Ca2 * g.BlockageFactor);
                double rm2_new = Math.Sqrt(A2 * (1.0 + ht) / (4.0 * Math.PI * (1.0 - ht)));

                if (Math.Abs(rm2_new - rm2) / Math.Max(rm2, 1e-9) < 1e-9)
                {
                    var s2 = new StationState();
                    FillStationGeometry(ref s2, rm2_new, ht, A2);

                    s2.U = U2;
                    s2.Ca = Ca2;
                    s2.Ct = Ct2;
                    s2.C = C2;
                    s2.Alpha = alpha2;

                    s2.Wa = Ca2;
                    s2.Wt = Wt2;
                    s2.W = W2;
                    s2.Beta = beta2;

                    s2.T = T2;
                    s2.P = P2;
                    s2.Rho = rho2;
                    s2.S = s2_static;
                    s2.H = h2;

                    s2.H0 = h02;
                    s2.T0 = T02;
                    s2.P0 = FlowProperties.ComputeTotalPressure(P2, T2, T02);

                    double T0rel2 = T2 + W2 * W2 / (2.0 * FlowProperties.Cp);
                    s2.T0Rel = T0rel2;
                    s2.P0Rel = FlowProperties.ComputeTotalPressure(P2, T2, T0rel2);

                    s2.Mach = FlowProperties.ComputeMach(C2, T2);
                    return s2;
                }

                rm2 = 0.5 * (rm2 + rm2_new);
            }

            // fallback minimal
            return new StationState { Radius = rm2 };
        }

        // =========================
        // STATION 3 : Stator Exit
        // =========================
        private StationState ComputeStation3(StageInput stg, MeanlineInputs g, StationState s2)
        {
            double ht = stg.HubTipRatio;
            double rm3 = s2.Radius;

            for (int iter = 0; iter < GeoMaxIter; iter++)
            {
                double U = s2.U;
                double Ca3 = stg.AxialVelocityRatio * s2.Ca;

                // Reaction constraint (deterministic):
                // R = 1 - (Ct2 + Ct3)/(2U) -> Ct3 = 2U(1-R) - Ct2
                double Ct3 = 2.0 * U * (1.0 - stg.Reaction) - s2.Ct;
                double alpha3 = Math.Atan2(Ct3, Ca3);
                double C3 = Math.Sqrt(Ca3 * Ca3 + Ct3 * Ct3);

                // h03 = h02
                double h03 = s2.H0;
                double T03 = h03 / FlowProperties.Cp;

                double h3 = h03 - 0.5 * C3 * C3;
                double T3 = h3 / FlowProperties.Cp;

                double Yst = EmpiricalModels.ComputeLossY_Stator(g.DiffusionFactor, g.ThicknessChordRatio);
                double dS23 = EmpiricalModels.ComputeEntropyRiseFromY(Yst);
                double s3_static = s2.S + dS23;

                double P3 = FlowProperties.ComputePressureFromEntropy(T3, s3_static);
                double rho3 = FlowProperties.ComputeDensity(P3, T3);

                double A3 = g.MassFlowRate / (rho3 * Ca3 * g.BlockageFactor);
                double rm3_new = Math.Sqrt(A3 * (1.0 + ht) / (4.0 * Math.PI * (1.0 - ht)));

                if (Math.Abs(rm3_new - rm3) / Math.Max(rm3, 1e-9) < 1e-9)
                {
                    var s3 = new StationState();
                    FillStationGeometry(ref s3, rm3_new, ht, A3);

                    s3.U = U;
                    s3.Ca = Ca3;
                    s3.Ct = Ct3;
                    s3.C = C3;
                    s3.Alpha = alpha3;

                    s3.T = T3;
                    s3.P = P3;
                    s3.Rho = rho3;
                    s3.S = s3_static;
                    s3.H = h3;

                    s3.H0 = h03;
                    s3.T0 = T03;
                    s3.P0 = FlowProperties.ComputeTotalPressure(P3, T3, T03);
                    s3.Mach = FlowProperties.ComputeMach(C3, T3);
                    return s3;
                }

                rm3 = 0.5 * (rm3 + rm3_new);
            }

            return new StationState { Radius = rm3 };
        }

        private BladeAngles ComputeBladeAngles(StageInput stg, MeanlineInputs g, StationState s1, StationState s2, StationState s3)
        {
            double sigma = EmpiricalModels.ComputeSolidityFromDiffusionFactor(g.DiffusionFactor);

            double i = EmpiricalModels.ComputeIncidence(g.DiffusionFactor, stg.Reaction, sigma);
            double d = EmpiricalModels.ComputeDeviation(g.DiffusionFactor, g.ThicknessChordRatio);

            double beta1_metal = s1.Beta - i;
            double beta2_metal = s2.Beta + d;

            double alpha2_metal = s2.Alpha - i;
            double alpha3_metal = s3.Alpha + d;

            double chordRotor = s1.Height / Math.Max(g.AspectRatio, 1e-9);
            double chordStator = s2.Height / Math.Max(g.AspectRatio, 1e-9);

            double camberRotor = beta1_metal - beta2_metal;
            double camberStator = alpha2_metal - alpha3_metal;

            return new BladeAngles
            {
                Solidity = sigma,
                Incidence = i,
                Deviation = d,
                RotorInletMetal = beta1_metal,
                RotorExitMetal = beta2_metal,
                StatorInletMetal = alpha2_metal,
                StatorExitMetal = alpha3_metal,
                ChordRotor = chordRotor,
                ChordStator = chordStator,
                CamberRotor = camberRotor,
                CamberStator = camberStator,
            };
        }

        private static void FillStationGeometry(ref StationState s, double rm, double hubTipRatio, double area)
        {
            double rt = 2.0 * rm / (1.0 + hubTipRatio);
            double rh = hubTipRatio * rt;

            s.Radius = rm;
            s.TipRadius = rt;
            s.HubRadius = rh;
            s.Area = area;
            s.Height = rt - rh;
        }

        private static void FillKinematicsAndThermo(ref StationState s, double omega, double Ca, double Ct, double alpha, double T0, double P0)
        {
            s.U = omega * s.Radius;

            s.Ca = Ca;
            s.Ct = Ct;
            s.C = Math.Sqrt(Ca * Ca + Ct * Ct);
            s.Alpha = alpha;

            s.Wa = Ca;
            s.Wt = s.U - Ct;
            s.W = Math.Sqrt(s.Wa * s.Wa + s.Wt * s.Wt);
            s.Beta = Math.Atan2(s.Wt, s.Wa);

            s.T0 = T0;
            s.P0 = P0;

            s.T = FlowProperties.ComputeStaticTemperature(T0, s.C);
            s.P = FlowProperties.ComputeStaticPressure(P0, s.T, T0);
            s.Rho = FlowProperties.ComputeDensity(s.P, s.T);

            s.S = FlowProperties.ComputeEntropy(s.T, s.P);
            s.H = FlowProperties.ComputeEnthalpy(s.T);
            s.H0 = FlowProperties.ComputeTotalEnthalpy(s.H, s.C);

            s.T0Rel = s.T + s.W * s.W / (2.0 * FlowProperties.Cp);
            s.P0Rel = FlowProperties.ComputeTotalPressure(s.P, s.T, s.T0Rel);

            s.Mach = FlowProperties.ComputeMach(s.C, s.T);
        }

        private sealed class BladeAngles
        {
            public double Solidity { get; set; }
            public double Incidence { get; set; }
            public double Deviation { get; set; }
            public double RotorInletMetal { get; set; }
            public double RotorExitMetal { get; set; }
            public double StatorInletMetal { get; set; }
            public double StatorExitMetal { get; set; }
            public double ChordRotor { get; set; }
            public double ChordStator { get; set; }
            public double CamberRotor { get; set; }
            public double CamberStator { get; set; }
        }
    }
}
