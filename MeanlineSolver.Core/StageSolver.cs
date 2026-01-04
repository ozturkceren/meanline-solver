using System;

namespace MeanlineSolver.Core
{
    public class StageSolver
    {
        /// <summary>
        /// Solve a single compressor stage at meanline level.
        /// For now, it computes Station 1 and Station 2 and basic stage performance.
        /// </summary>
        public StageResult SolveSingleStage(StageInput stage, MeanlineInputs global)
        {
            // Station 1: rotor inlet (uses global inlet conditions for first stage)
            var st1 = ComputeStation1(stage, global);

            // Station 2: rotor exit
            var st2 = ComputeStation2(stage, global, st1);

            var st3 = ComputeStation3(stage, global, st2);

            // Stage-level performance
            double deltaH0 = st2.H0 - st1.H0;                                 // total enthalpy rise
            double phi = VelocityTriangle.ComputePhi(st1.Ca, st1.U);         // flow coefficient at inlet
            double psi = VelocityTriangle.ComputePsi(deltaH0, st1.U);        // loading coefficient

            // Isentropic pressure ratio (ideal)
            double prIs = FlowProperties.ComputeIsentropicPressureRatio(deltaH0, st1.T0);

            // Loss-based efficiency estimate
            double totalLoss = EmpiricalModels.ComputeTotalLoss(
                stage.ThicknessChordRatio,
                stage.Solidity,
                stage.HubTipRatio,
                global.TipClearance,
                st1.Height);

            double etaStage = EmpiricalModels.EstimateEfficiencyFromLoss(totalLoss);

            // Actual stage pressure ratio (very simple mapping)
            double prActual = 1.0 + etaStage * (prIs - 1.0);

            // Placeholder blade angles (ileride ComputeBladeAngles ile dolduracağız)
            var blades = ComputeBladeAngles(stage, st1, st2, st3);

            var result = new StageResult
            {
                Station1 = st1,
                Station2 = st2,
                Station3 = st3,

                DeltaH0 = deltaH0,
                PressureRatio = prActual,
                Efficiency = etaStage,

                FlowCoefficient = phi,
                LoadingCoefficient = psi,
                Reaction = stage.Reaction,

                RotorInletAngle = blades.rotorInletBladeAngle,
                RotorExitAngle = blades.rotorExitBladeAngle,
                StatorInletAngle = blades.statorInletBladeAngle,
                StatorExitAngle = blades.statorExitBladeAngle,

                Incidence = blades.incidence,
                Deviation = blades.deviation,
                Solidity = blades.solidity
            };

            // CHORD LENGTHS
            result.ChordLengthRotor = st1.Height / blades.solidity;
            result.ChordLengthStator = st2.Height / blades.solidity;

            // CAMBER ANGLES
            result.CamberAngleRotor = blades.rotorInletBladeAngle - blades.rotorExitBladeAngle;
            result.CamberAngleStator = blades.statorInletBladeAngle - blades.statorExitBladeAngle;

            return result;

        }

        // -------------------------------------------------
        //  STATION 1 : ROTOR INLET  (2.1.5.1)
        // -------------------------------------------------
        public StationState ComputeStation1(StageInput stage, MeanlineInputs global)
        {
            StationState st = new StationState();

            // 1) Mean radius
            // HubTipRatio = rh / rt  (global'den alıyoruz)
            double rhFrac = global.InletHubTipRatio; // 0.8 gibi bir oran
            double rt = global.TipRadius;                         // ilk stage için normalize tip radius
            double rh = rhFrac * rt;
            double rm = 0.5 * (rh + rt);

            st.Radius = rm;

            // Blade speed: U = ω r  (ω: rad/s)
            double omega = global.RotationalSpeed * 2.0 * Math.PI / 60.0;
            st.U = omega * rm;

            // 2) Axial velocity: Ca1 = φ * U
            st.Ca = stage.FlowCoefficient * st.U;

            // 3) Tangential velocity: inlet swirl ~ 0 varsayımı
            st.Ct = 0.0;

            // 4) Absolute velocity magnitude
            st.C = VelocityTriangle.ComputeAbsoluteVelocity(st.Ca, st.Ct);

            // 5) Flow angle α1
            st.Alpha = VelocityTriangle.ComputeAlpha(st.Ca, st.Ct);

            // 6) Total → static: T1, P1
            st.T0 = global.InletTotalTemperature;
            st.P0 = global.InletTotalPressure;

            st.T = FlowProperties.ComputeStaticTemperature(st.T0, st.C);
            st.P = FlowProperties.ComputeStaticPressure(st.P0, st.T, st.T0);

            // 7) Density, area, height
            st.Rho = FlowProperties.ComputeDensity(st.P, st.T);

            st.Area = global.MassFlowRate / (st.Rho * st.Ca);
            st.Height = st.Area / (2.0 * Math.PI * st.Radius);

            // 8) Relative velocity & β1
            (st.Wa, st.Wt, st.W) =
                VelocityTriangle.ComputeRelativeVelocity(st.Ca, st.Ct, st.U);

            st.Beta = VelocityTriangle.ComputeBeta(st.Wa, st.Wt);

            // 9) Enthalpy, total enthalpy, Mach
            st.H = FlowProperties.ComputeEnthalpy(st.T);
            st.H0 = FlowProperties.ComputeTotalEnthalpy(st.H, st.C);
            st.Mach = FlowProperties.ComputeMach(st.C, st.T);

            return st;
        }

        // -------------------------------------------------
        //  STATION 2 : ROTOR EXIT  (basitleştirilmiş 2.1.5.2)
        // -------------------------------------------------
        public StationState ComputeStation2(StageInput stage, MeanlineInputs global, StationState st1)
        {
            StationState st2 = new StationState();

            // 1) Radius and blade speed (şimdilik rm2 = rm1 alıyoruz)
            st2.Radius = st1.Radius;
            double omega = global.RotationalSpeed * 2.0 * Math.PI / 60.0;
            st2.U = omega * st2.Radius;

            // 2) Axial velocity: Ca2 = Ca1 * AVR_rotor
            st2.Ca = VelocityTriangle.ComputeNextAxialVelocity(
                st1.Ca,
                stage.AxialVelocityRatioRotor);

            // 3) Tangential velocity: Ct2 from loading coefficient ψ
            // ψ = Δh0 / U^2 = (Ct2 - Ct1) / U  → Ct2 = Ct1 + ψ U
            st2.Ct = VelocityTriangle.ComputeCt2(st1.Ct, stage.LoadingCoefficient, st2.U);

            // 4) Absolute velocity C2 and flow angle α2
            st2.C = VelocityTriangle.ComputeAbsoluteVelocity(st2.Ca, st2.Ct);
            st2.Alpha = VelocityTriangle.ComputeAlpha(st2.Ca, st2.Ct);

            // 5) Rotor work and total enthalpy
            double deltaH0 = VelocityTriangle.ComputeDeltaH0(st2.U, st2.Ct, st1.Ct);
            // h0,2 = h0,1 + Δh0
            st2.H0 = st1.H0 + deltaH0;

            // Corresponding total temperature T0,2 from h0 = Cp T0
            st2.T0 = st2.H0 / FlowProperties.Cp;

            // 6) Static temperature (T2) and pressure (P2)
            st2.T = FlowProperties.ComputeStaticTemperature(st2.T0, st2.C);

            // Isentropic pressure ratio for rotor (ideal)
            double prIs = FlowProperties.ComputeIsentropicPressureRatio(deltaH0, st1.T0);

            // Loss-based efficiency estimate for rotor part
            double totalLoss = EmpiricalModels.ComputeTotalLoss(
                stage.ThicknessChordRatio,
                stage.Solidity,
                stage.HubTipRatio,
                global.TipClearance,
                st1.Height); // clearance loss ~ inlet height

            double etaRotor = EmpiricalModels.EstimateEfficiencyFromLoss(totalLoss);

            // Actual total pressure ratio across rotor
            double prActual = 1.0 + etaRotor * (prIs - 1.0);

            st2.P0 = st1.P0 * prActual;

            // Static pressure from total pressure
            st2.P = FlowProperties.ComputeStaticPressure(st2.P0, st2.T, st2.T0);

            // 7) Density, area, height
            st2.Rho = FlowProperties.ComputeDensity(st2.P, st2.T);

            st2.Area = global.MassFlowRate / (st2.Rho * st2.Ca);
            st2.Height = st2.Area / (2.0 * Math.PI * st2.Radius);

            // 8) Relative velocity and β2
            (st2.Wa, st2.Wt, st2.W) =
                VelocityTriangle.ComputeRelativeVelocity(st2.Ca, st2.Ct, st2.U);

            st2.Beta = VelocityTriangle.ComputeBeta(st2.Wa, st2.Wt);

            // 9) Mach number and static enthalpy
            st2.H = FlowProperties.ComputeEnthalpy(st2.T);
            st2.Mach = FlowProperties.ComputeMach(st2.C, st2.T);

            return st2;
        }
        public StationState ComputeStation3(StageInput stage, MeanlineInputs global, StationState st2)
        {
            StationState st3 = new StationState();

            // 1) Radius and blade speed (şimdilik sabit mean radius)
            st3.Radius = st2.Radius;
            double omega = global.RotationalSpeed * 2.0 * Math.PI / 60.0;
            st3.U = omega * st3.Radius;

            // 2) Axial velocity: Ca3 = Ca2 * AVR_stator
            st3.Ca = VelocityTriangle.ComputeNextAxialVelocity(
                st2.Ca,
                stage.AxialVelocityRatioStator);

            // 3) Flow angle α3 (simplified model)
            double alpha3 = stage.InletFlowAngle * Math.PI / 180.0; // deg → rad
            st3.Alpha = alpha3;

            // 4) Tangential Ct3 = Ca3 * tan(alpha3)
            st3.Ct = st3.Ca * Math.Tan(alpha3);

            // 5) Absolute velocity magnitude
            st3.C = VelocityTriangle.ComputeAbsoluteVelocity(st3.Ca, st3.Ct);

            // 6) Total enthalpy conserved
            st3.H0 = st2.H0;
            st3.T0 = st2.T0;
            st3.P0 = st2.P0;

            // 7) Static properties
            st3.T = FlowProperties.ComputeStaticTemperature(st3.T0, st3.C);
            st3.P = FlowProperties.ComputeStaticPressure(st3.P0, st3.T, st3.T0);
            st3.Rho = FlowProperties.ComputeDensity(st3.P, st3.T);

            // 8) Area and blade height
            st3.Area = global.MassFlowRate / (st3.Rho * st3.Ca);
            st3.Height = st3.Area / (2.0 * Math.PI * st3.Radius);

            // 9) Relative velocity and β3
            (st3.Wa, st3.Wt, st3.W) =
                VelocityTriangle.ComputeRelativeVelocity(st3.Ca, st3.Ct, st3.U);

            st3.Beta = VelocityTriangle.ComputeBeta(st3.Wa, st3.Wt);

            // 10) Enthalpy & Mach
            st3.H = FlowProperties.ComputeEnthalpy(st3.T);
            st3.Mach = FlowProperties.ComputeMach(st3.C, st3.T);

            return st3;
        }
        public (double rotorInletBladeAngle,
        double rotorExitBladeAngle,
        double statorInletBladeAngle,
        double statorExitBladeAngle,
        double incidence,
        double deviation,
        double solidity)
    ComputeBladeAngles(StageInput stage, StationState st1, StationState st2, StationState st3)
        {
            // Başlangıç solidity tahmini
            double sigma = stage.Solidity;
            double oldSigma;

            double incidence = 0.0;
            double deviation = 0.0;

            // Maks 20 iterasyon
            for (int iter = 0; iter < 20; iter++)
            {
                oldSigma = sigma;

                // 1) IDEAL açılar (velocity triangles)
                double alpha1 = st1.Alpha;
                double beta1 = st1.Beta;
                double alpha2 = st2.Alpha;
                double beta2 = st2.Beta;
                double alpha3 = st3.Alpha;

                // 2) INCIDENCE
                incidence = EmpiricalModels.ComputeIncidence(
                    stage.FlowCoefficient,
                    stage.Reaction,
                    sigma);

                // 3) DEVIATION
                deviation = EmpiricalModels.ComputeDeviation(
                    stage.ThicknessChordRatio,
                    sigma);

                // 4) BLADE AÇILARI
                double beta1Blade = beta1 - incidence;     // rotor inlet blade angle
                double beta2Blade = beta2 + deviation;     // rotor exit blade angle

                double alpha2Blade = alpha2;               // stator inlet
                double alpha3Blade = alpha3 + deviation;   // stator exit

                // 5) SOLIDITY UPDATE  (çok basit model)
                // sigma_new = sigma_old * (1 + deviation * 0.5)
                sigma = sigma * (1.0 + 0.5 * Math.Abs(deviation));

                // KONVERJANS kontrolü
                if (Math.Abs(sigma - oldSigma) < 1e-5)
                {
                    return (beta1Blade, beta2Blade, alpha2Blade, alpha3Blade,
                            incidence, deviation, sigma);
                }
            }

            // Konverjans sağlanmasa bile en son değerleri döndür
            return (st1.Beta, st2.Beta, st2.Alpha, st3.Alpha, incidence, deviation, sigma);
        }

    }
}
