using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace WindowsFormsApplication1
{
    
    class State
    {
        //state contains a reduced matrix, a lowerbound, a route traveled and the list of cities visited
        public double[,] reducedMatrix;
        public double lowerBound;
        public List<int> route;
        public bool[] citiesVisited;

        //constructs a state
        public State(double[,] reducedMatrix, double lowerBound, List<int> route, bool[] citiesVisited)
        {
            this.reducedMatrix = reducedMatrix;
            this.lowerBound = lowerBound;
            this.route = route;
            this.citiesVisited = citiesVisited;
        }
        
        //returns the priority value to the priority queue
        //this favors going deep in the tree as 1000 pixels traveled to encourage depth over branching
        public double getPriorityValue()
        {
            double val = lowerBound - route.Count * 1000;
            return lowerBound - route.Count * 1000;
        }
    }
}
