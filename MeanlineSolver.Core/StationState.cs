using System;

namespace MeanlineSolver.Core
{
    public class StationState
    {
        // Velocities
        public double Ca { get; set; }
        public double Ct { get; set; }
        public double C { get; set; }

        public double Wa { get; set; }
        public double Wt { get; set; }
        public double W { get; set; }

        // Angles
        public double Alpha { get; set; }
        public double Beta { get; set; }

        // Thermodynamic properties
        public double T { get; set; }
        public double P { get; set; }
        public double Rho { get; set; }
        public double H { get; set; }
        public double S { get; set; }

        // Relative totals (rotor completeness)
        public double T0Rel { get; set; }
        public double P0Rel { get; set; }

        // Absolute totals
        public double T0 { get; set; }
        public double P0 { get; set; }
        public double H0 { get; set; }

        // Geometry
        public double Radius { get; set; }     // rm
        public double HubRadius { get; set; }  // rh
        public double TipRadius { get; set; }  // rt
        public double Height { get; set; }     // rt - rh
        public double Area { get; set; }       // flow area

        public double Mach { get; set; }
        public double U { get; set; }
    }
}
