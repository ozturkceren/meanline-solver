using System;

namespace MeanlineSolver.Core
{
    public class StageSolver
    {
        // Default numerical settings (can be lifted to MeanlineInputs later)
        private const int OuterMaxIter = 40;
        private const int GeoMaxIter = 25;
        private const double OuterTolReaction = 1e-4;
        private const double OuterTolPr = 5e-4;
        private const double GeoTolRm = 1e-6;

        /// <summary>
        /// Solve a single compressor stage at meanline level.
        ///
        /// This implementation follows your notebook structure at a practical level:
        /// - outer loop: adjust ψ (loading) and α3 (stator exit angle) to match target reaction and target PR/stage
        /// - inner loop (Station 2/3): iterate mean radius via continuity (annulus geometry) until rm converges
        ///
        /// IMPORTANT: This is still a compact meanline model; your EmpiricalModels are placeholders.
        /// </summary>
        public StageResult SolveSingleStage(StageInput stage, MeanlineInputs global)
        {
            if (stage == null) throw new ArgumentNullException(nameof(stage));
            if (global == null) throw new ArgumentNullException(nameof(global));

            // Stage targets
            int nStages = Math.Max(global.Stages?.Count ?? 1, 1);
            double targetPrStage = Math.Pow(Math.Max(global.PressureRatio, 1.0), 1.0 / nStages);
            double targetReaction = stage.Reaction;

            // Design variables we will adjust
            double psi = stage.LoadingCoefficient;
            double alpha3Deg = stage.StatorExitFlowAngle;
            if (Math.Abs(alpha3Deg) < 1e-12) alpha3Deg = 0.0; // explicit

            StationState st1 = null!;
            StationState st2 = null!;
            StationState st3 = null!;

            double prActual = 1.0;
            double etaStage = 0.9;
            double deltaH0 = 0.0;
            double reactionCalc = stage.Reaction;

            // Outer iteration
            for (int it = 0; it < OuterMaxIter; it++)
            {
                // --- Station 1 (rotor inlet) ---
                st1 = ComputeStation1(stage, global);

                // --- Station 2 (rotor exit) with geometry iteration ---
                st2 = ComputeStation2(stage, global, st1, psi);

                // --- Station 3 (stator exit) with geometry iteration ---
                st3 = ComputeStation3(stage, global, st2, alpha3Deg);

                // Stage total enthalpy rise is rotor work in this compact model
                deltaH0 = st2.H0 - st1.H0;

                // Loss / efficiency
                double totalLoss = EmpiricalModels.ComputeTotalLoss(
                    stage.ThicknessChordRatio,
                    stage.Solidity,
                    stage.HubTipRatio,
                    global.TipClearance,
                    Math.Max(st1.Height, 1e-9));

                etaStage = EmpiricalModels.EstimateEfficiencyFromLoss(totalLoss);

                // Map to total pressure ratio (stage)
                double prIs = FlowProperties.ComputeIsentropicPressureRatio(deltaH0, st1.T0);
                prActual = 1.0 + etaStage * (prIs - 1.0);

                // Reaction estimate (from static enthalpy split)
                reactionCalc = ComputeReactionFromStates(st1, st2, st3);

                // Check convergence
                double errR = targetReaction - reactionCalc;
                double errPr = targetPrStage - prActual;
                if (Math.Abs(errR) < OuterTolReaction && Math.Abs(errPr) < OuterTolPr)
                    break;

                // --- Update ψ to hit PR target (dominant driver) ---
                // Required Δh0 for target PR (inverting isentropic relation with efficiency)
                double term = Math.Pow(Math.Max(targetPrStage, 1.0), (FlowProperties.Gamma - 1.0) / FlowProperties.Gamma);
                double deltaH0ReqIs = FlowProperties.Cp * st1.T0 * (term - 1.0);
                double deltaH0Req = deltaH0ReqIs / Math.Max(etaStage, 1e-6);
                double psiReq = deltaH0Req / (st1.U * st1.U);

                // Relaxation
                psi = 0.7 * psi + 0.3 * psiReq;
                psi = Clamp(psi, 0.05, 1.2);

                // --- Update α3 to hit reaction target (numerical derivative) ---
                // Use a small perturbation around current α3
                double dAlpha = 0.5; // degrees
                double r0 = reactionCalc;
                double r1 = ComputeReactionWithAlpha3(stage, global, st1, psi, alpha3Deg + dAlpha);
                double dRdA = (r1 - r0) / (dAlpha * Math.PI / 180.0); // per rad

                if (Math.Abs(dRdA) > 1e-8)
                {
                    double stepRad = errR / dRdA;
                    double stepDeg = stepRad * 180.0 / Math.PI;
                    alpha3Deg = alpha3Deg + Clamp(stepDeg, -2.0, 2.0); // limit step size
                }
                else
                {
                    // Fallback: gentle bias
                    alpha3Deg += Clamp(errR * 2.0, -1.0, 1.0);
                }

                alpha3Deg = Clamp(alpha3Deg, -25.0, 25.0);
            }

            // Blade angles (iterative incidence/deviation + solidity update)
            var blades = IterateAngles.Solve(
                phi: stage.FlowCoefficient,
                reactionTarget: targetReaction,
                thicknessChordRatio: stage.ThicknessChordRatio,
                solidityInitial: stage.Solidity,
                alpha1: st1.Alpha,
                beta1: st1.Beta,
                alpha2: st2.Alpha,
                beta2: st2.Beta,
                alpha3: st3.Alpha);

            var result = new StageResult
            {
                Station1 = st1,
                Station2 = st2,
                Station3 = st3,

                DeltaH0 = deltaH0,
                PressureRatio = prActual,
                Efficiency = etaStage,

                FlowCoefficient = stage.FlowCoefficient,
                LoadingCoefficient = psi,
                Reaction = reactionCalc,

                RotorInletAngle = blades.RotorInletMetal,
                RotorExitAngle = blades.RotorExitMetal,
                StatorInletAngle = blades.StatorInletMetal,
                StatorExitAngle = blades.StatorExitMetal,

                Incidence = blades.Incidence,
                Deviation = blades.Deviation,
                Solidity = blades.Solidity
            };

            // CHORD LENGTHS (very basic: c = h/σ)
            result.ChordLengthRotor = st1.Height / Math.Max(blades.Solidity, 1e-9);
            result.ChordLengthStator = st2.Height / Math.Max(blades.Solidity, 1e-9);

            // CAMBER ANGLES
            result.CamberAngleRotor = blades.RotorInletMetal - blades.RotorExitMetal;
            result.CamberAngleStator = blades.StatorInletMetal - blades.StatorExitMetal;

            return result;
        }

        // -------------------------------------------------
        //  STATION 1 : ROTOR INLET
        // -------------------------------------------------
        public StationState ComputeStation1(StageInput stage, MeanlineInputs global)
        {
            var st = new StationState();

            // Mean radius from inlet tip radius and inlet hub-tip ratio
            double rt = global.TipRadius;
            double rh = global.InletHubTipRatio * rt;
            double rm = 0.5 * (rh + rt);
            st.Radius = rm;

            // Blade speed
            double omega = global.RotationalSpeed * 2.0 * Math.PI / 60.0;
            st.U = omega * rm;

            // Axial velocity from φ
            st.Ca = stage.FlowCoefficient * st.U;

            // Optional inlet swirl from α1 (deg)
            double alpha1 = stage.InletFlowAngle * Math.PI / 180.0;
            st.Alpha = alpha1;
            st.Ct = st.Ca * Math.Tan(alpha1);

            st.C = VelocityTriangle.ComputeAbsoluteVelocity(st.Ca, st.Ct);

            // Total conditions
            st.T0 = global.InletTotalTemperature;
            st.P0 = global.InletTotalPressure;

            // Static
            st.T = FlowProperties.ComputeStaticTemperature(st.T0, st.C);
            st.P = FlowProperties.ComputeStaticPressure(st.P0, st.T, st.T0);
            st.Rho = FlowProperties.ComputeDensity(st.P, st.T);

            // Continuity area (include blockage)
            double B = (global.BlockageFactor <= 0.0) ? 1.0 : global.BlockageFactor;
            st.Area = global.MassFlowRate / (st.Rho * st.Ca * B);
            st.Height = st.Area / (2.0 * Math.PI * st.Radius);

            // Relative
            (st.Wa, st.Wt, st.W) = VelocityTriangle.ComputeRelativeVelocity(st.Ca, st.Ct, st.U);
            st.Beta = VelocityTriangle.ComputeBeta(st.Wa, st.Wt);

            // Energy
            st.H = FlowProperties.ComputeEnthalpy(st.T);
            st.H0 = FlowProperties.ComputeTotalEnthalpy(st.H, st.C);
            st.Mach = FlowProperties.ComputeMach(st.C, st.T);

            return st;
        }

        // -------------------------------------------------
        //  STATION 2 : ROTOR EXIT (with rm iteration)
        // -------------------------------------------------
        public StationState ComputeStation2(StageInput stage, MeanlineInputs global, StationState st1, double psi)
        {
            var st2 = new StationState();

            double rm = st1.Radius;
            double omega = global.RotationalSpeed * 2.0 * Math.PI / 60.0;

            // Iteration on mean radius via area relation (annulus with given hub-tip ratio)
            for (int k = 0; k < GeoMaxIter; k++)
            {
                double rmOld = rm;

                st2.Radius = rm;
                st2.U = omega * rm;

                // Axial velocity
                st2.Ca = VelocityTriangle.ComputeNextAxialVelocity(st1.Ca, stage.AxialVelocityRatioRotor);

                // Tangential velocity (notebook form)
                st2.Ct = VelocityTriangle.ComputeCt2(st1.Ct, psi, st2.U, st1.Radius, st2.Radius);

                // Absolute
                st2.C = VelocityTriangle.ComputeAbsoluteVelocity(st2.Ca, st2.Ct);
                st2.Alpha = VelocityTriangle.ComputeAlpha(st2.Ca, st2.Ct);

                // Euler work
                double deltaH0 = VelocityTriangle.ComputeDeltaH0(st1.U, st2.U, st1.Ct, st2.Ct);
                st2.H0 = st1.H0 + deltaH0;
                st2.T0 = st2.H0 / FlowProperties.Cp;

                // Loss-based efficiency (rotor part approximation)
                double totalLoss = EmpiricalModels.ComputeTotalLoss(
                    stage.ThicknessChordRatio,
                    stage.Solidity,
                    stage.HubTipRatio,
                    global.TipClearance,
                    Math.Max(st1.Height, 1e-9));
                double etaRotor = EmpiricalModels.EstimateEfficiencyFromLoss(totalLoss);

                // Total pressure rise (mapping)
                double prIs = FlowProperties.ComputeIsentropicPressureRatio(deltaH0, st1.T0);
                double prActual = 1.0 + etaRotor * (prIs - 1.0);
                st2.P0 = st1.P0 * prActual;

                // Static
                st2.T = FlowProperties.ComputeStaticTemperature(st2.T0, st2.C);
                st2.P = FlowProperties.ComputeStaticPressure(st2.P0, st2.T, st2.T0);
                st2.Rho = FlowProperties.ComputeDensity(st2.P, st2.T);

                // Continuity area (include blockage)
                double B = (global.BlockageFactor <= 0.0) ? 1.0 : global.BlockageFactor;
                st2.Area = global.MassFlowRate / (st2.Rho * st2.Ca * B);

                // Update rm from annulus relation (hub-tip ratio fixed for this stage)
                rm = ComputeMeanRadiusFromAreaAndHubTip(st2.Area, stage.HubTipRatio);
                if (Math.Abs(rm - rmOld) < GeoTolRm)
                {
                    break;
                }
            }

            st2.Height = st2.Area / (2.0 * Math.PI * st2.Radius);

            // Relative
            (st2.Wa, st2.Wt, st2.W) = VelocityTriangle.ComputeRelativeVelocity(st2.Ca, st2.Ct, st2.U);
            st2.Beta = VelocityTriangle.ComputeBeta(st2.Wa, st2.Wt);

            // Static enthalpy & Mach
            st2.H = FlowProperties.ComputeEnthalpy(st2.T);
            st2.Mach = FlowProperties.ComputeMach(st2.C, st2.T);

            return st2;
        }

        // -------------------------------------------------
        //  STATION 3 : STATOR EXIT (with rm iteration)
        // -------------------------------------------------
        public StationState ComputeStation3(StageInput stage, MeanlineInputs global, StationState st2, double alpha3Deg)
        {
            var st3 = new StationState();

            double rm = st2.Radius;
            double omega = global.RotationalSpeed * 2.0 * Math.PI / 60.0;
            double alpha3 = alpha3Deg * Math.PI / 180.0;

            for (int k = 0; k < GeoMaxIter; k++)
            {
                double rmOld = rm;

                st3.Radius = rm;
                st3.U = omega * rm;

                // Axial velocity
                st3.Ca = VelocityTriangle.ComputeNextAxialVelocity(st2.Ca, stage.AxialVelocityRatioStator);

                // Set stator exit angle
                st3.Alpha = alpha3;
                st3.Ct = st3.Ca * Math.Tan(alpha3);
                st3.C = VelocityTriangle.ComputeAbsoluteVelocity(st3.Ca, st3.Ct);

                // Total enthalpy approximately conserved across stator in meanline
                st3.H0 = st2.H0;
                st3.T0 = st2.T0;
                st3.P0 = st2.P0;

                // Static
                st3.T = FlowProperties.ComputeStaticTemperature(st3.T0, st3.C);
                st3.P = FlowProperties.ComputeStaticPressure(st3.P0, st3.T, st3.T0);
                st3.Rho = FlowProperties.ComputeDensity(st3.P, st3.T);

                // Continuity area (include blockage)
                double B = (global.BlockageFactor <= 0.0) ? 1.0 : global.BlockageFactor;
                st3.Area = global.MassFlowRate / (st3.Rho * st3.Ca * B);

                // Update rm (same hub-tip ratio assumption)
                rm = ComputeMeanRadiusFromAreaAndHubTip(st3.Area, stage.HubTipRatio);
                if (Math.Abs(rm - rmOld) < GeoTolRm)
                    break;
            }

            st3.Height = st3.Area / (2.0 * Math.PI * st3.Radius);

            // Relative (not used for compressor stator, but kept for completeness)
            (st3.Wa, st3.Wt, st3.W) = VelocityTriangle.ComputeRelativeVelocity(st3.Ca, st3.Ct, st3.U);
            st3.Beta = VelocityTriangle.ComputeBeta(st3.Wa, st3.Wt);

            st3.H = FlowProperties.ComputeEnthalpy(st3.T);
            st3.Mach = FlowProperties.ComputeMach(st3.C, st3.T);

            return st3;
        }

        // -------------------------------------------------
        //  Helpers
        // -------------------------------------------------

        private static double ComputeMeanRadiusFromAreaAndHubTip(double area, double hubTipRatio)
        {
            // Derivation (hub-tip ratio = rh/rt, rm=(rh+rt)/2) gives:
            // A = 4π rm^2 * (1 - ht) / (1 + ht)
            double ht = Clamp(hubTipRatio, 0.2, 0.99);
            double denom = 4.0 * Math.PI * (1.0 - ht);
            double num = area * (1.0 + ht);
            if (num <= 0.0 || denom <= 0.0)
                return 0.0;
            return Math.Sqrt(num / denom);
        }

        private static double ComputeReactionFromStates(StationState st1, StationState st2, StationState st3)
        {
            // R = (h2 - h1) / (h3 - h1)
            double denom = (st3.H - st1.H);
            if (Math.Abs(denom) < 1e-12)
                return 0.5;
            return (st2.H - st1.H) / denom;
        }

        private double ComputeReactionWithAlpha3(StageInput stage, MeanlineInputs global, StationState st1, double psi, double alpha3Deg)
        {
            var st2 = ComputeStation2(stage, global, st1, psi);
            var st3 = ComputeStation3(stage, global, st2, alpha3Deg);
            return ComputeReactionFromStates(st1, st2, st3);
        }

        private static double Clamp(double x, double lo, double hi)
        {
            if (x < lo) return lo;
            if (x > hi) return hi;
            return x;
        }
    }
}
