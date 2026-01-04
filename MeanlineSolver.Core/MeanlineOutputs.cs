using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeanlineSolver.Core
{
    public class MeanlineOutputs
    {
        public double[] FlowInletAngles { get; set; }
        public double[] FlowOutletAngles { get; set; }
        public double[] BladeInletAngles { get; set; }
        public double[] BladeOutletAngles { get; set; }
        public double[] FlowInletVelocities { get; set; }
        public double[] FlowOutletVelocities { get; set; }
        public double[] IncidenceAngles { get; set; }
        public double[] DeviationAngles { get; set; }
        public double[] CamberAngles { get; set; }
        public double[] PitchAngles { get; set; }
        public double[] ChordLengths { get; set; }
        public double[] Solidities { get; set; }
        public double[] HubRadii { get; set; }
        public double[] TipRadii { get; set; }
    }
}
