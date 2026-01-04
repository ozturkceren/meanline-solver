using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MeanlineSolver.Core
{
    public class Stage
    {
        public int Index { get; set; }
        public double Reaction { get; set; }
        public double LoadCoefficient { get; set; }
        public double Solidity { get; set; }
        public double Chord { get; set; }
    }

}
