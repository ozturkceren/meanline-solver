using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeanlineSolver.Core
{
    public class StageInput
    {
        // Fundamental meanline parameters
        public double Reaction { get; set; }            // Λ
        public double FlowCoefficient { get; set; }     // φ
        public double LoadingCoefficient { get; set; }  // ψ

        // Geometry
        public double HubTipRatio { get; set; }         // r_h / r_t

        // Axial velocity evolution (pdf: AVR rotor & stator)
        public double AxialVelocityRatioRotor { get; set; }   // Ca2 = Ca1 * AVR_rotor
        public double AxialVelocityRatioStator { get; set; }  // Ca3 = Ca2 * AVR_stator

        // Blade geometry parameters
        public double Solidity { get; set; }            // σ = c/s
        public double ThicknessChordRatio { get; set; } // t/c

        // Inlet/exit flow conditions (optional swirl)

        /// <summary>
        /// Inlet absolute flow angle α1 in degrees (positive in direction of rotor rotation).
        /// If omitted in JSON, it defaults to 0.
        /// </summary>
        public double InletFlowAngle { get; set; }      // α1 [deg]

        /// <summary>
        /// Stator exit absolute flow angle α3 in degrees.
        /// If omitted in JSON, it defaults to 0.
        /// </summary>
        public double StatorExitFlowAngle { get; set; } // α3 [deg]
    }
}
