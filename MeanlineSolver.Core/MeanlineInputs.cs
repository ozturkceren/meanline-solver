using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeanlineSolver.Core
{
    public class MeanlineInputs
    {
        public string CompressorType { get; set; }
        public double MassFlowRate { get; set; }
        public int NumberOfStages { get; set; }
        public double PressureRatio { get; set; }
        public double RotationalSpeed { get; set; }
        public double AspectRatio { get; set; }
        public double[] StageReactions { get; set; }
        public double[] StageLoadDistribution { get; set; }
        public double TipClearance { get; set; }
        public double TipRadius { get; set; }
        public double ThicknessChordRatio { get; set; } // t/c
        public double BlockageFactor { get; set; }
        public double DiffusionFactor { get; set; }
        public double AxialVelocityRatio { get; set; }
        public double InletFlowAngle { get; set; }
        public double InletFlowCoefficient { get; set; }
        public double InletHubTipRatio { get; set; }
        public double InletTotalTemperature { get; set; }
        public double InletTotalPressure { get; set; }

        public List<StageInput> Stages { get; set; }
    }
}
