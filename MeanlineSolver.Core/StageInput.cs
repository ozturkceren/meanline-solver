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

        // Inlet flow conditions (optional swirl)
        public double InletFlowAngle { get; set; }      // α1, default = 0
    }
}
