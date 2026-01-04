 using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeanlineSolver.Core
{
        public class StationState
        {
            // Velocities
            public double Ca { get; set; }        // Axial velocity
            public double Ct { get; set; }        // Tangential velocity
            public double C { get; set; }         // Absolute velocity magnitude

            public double Wa { get; set; }        // Relative axial velocity
            public double Wt { get; set; }        // Relative tangential velocity
            public double W { get; set; }         // Relative velocity magnitude

            // Angles
            public double Alpha { get; set; }     // Flow angle
            public double Beta { get; set; }      // Relative flow angle

            // Thermodynamic properties
            public double T { get; set; }         // Static temperature
            public double P { get; set; }         // Static pressure
            public double Rho { get; set; }       // Density
            public double H { get; set; }         // Static enthalpy
            public double S { get; set; }         // Static entropy

            // Total (stagnation) properties
            public double T0 { get; set; }
            public double P0 { get; set; }
            public double H0 { get; set; }

            // Geometry
            public double Radius { get; set; }    // Mean radius
            public double Height { get; set; }    // Blade height
            public double Area { get; set; }      // Flow area

            // Derived quantities
            public double Mach { get; set; }
            public double U { get; set; }         // Blade speed at station
        }
    }
