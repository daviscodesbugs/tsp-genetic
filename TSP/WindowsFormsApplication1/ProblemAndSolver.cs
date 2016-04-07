using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;



namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                {
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
                    Cities[i].position = i;
                }
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                {
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
                    Cities[i].position = i;
                }
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            //the results that will be returned
            string[] results = new string[3];

            
            Stopwatch timer = new Stopwatch();
            TimeSpan endTime = new TimeSpan(0, 0, time_limit/1000);
            int numSolutionsFound = 0;
            timer.Start();

            //find our bssf and start at node 0
            //Uses a greedy algorithm that takes the closest node and continues on.
            //If in the hard version we end up getting stuck with an infinite route, we start on a different node
            //and run the greedy algorithm again.
            //If none of the nodes produces a viable solution, we run the default algoritm which finds a random route.
            List<int> route;
            //min keeps track of the smallest route we've found so far
            double min;
            //already gone keeps track of what cities we have already gone to to not look at those again
            bool[] citiesVisited;
            double tempCost;
            int tempIndex;
            int currentCity = -1;
            //if fail is true then we found an infinite route and nee
            bool fail;
            do //it is unlikely but we could have to run through all cities, and then go to the default solution so this could take speed O(n^3) with space complexity of O(n) for arrays and list
            {

                fail = false;
                currentCity++;
                citiesVisited = new bool[Cities.Length];
                citiesVisited[currentCity] = true;
                route = new List<int>();
                route.Add(currentCity);

                if (currentCity >= Cities.Length) //if all greedy attempts don't solve it, use the default solution
                {
                    defaultSolveProblem();
                    break;
                }
                //start at city 0 and find a route to all cities
                //we connect to every city and run through this n times which contains a O(n) method so it is O(n^2) with the space complexity of O(n) for arrays and list
                for (int i = 0; i < Cities.Length - 1; i++)
                {
                    min = double.PositiveInfinity;
                    tempIndex = -1;
                    //look at the route from our currentCity all available cities
                    //on average, we will run through this for loop half the time or O(n/2) which is O(n) speed complexity with the same size complexity (for cities visited)
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        //skip the routes where we have already gone
                        if (citiesVisited[j])
                            continue;
                        tempCost = Cities[currentCity].costToGetTo(Cities[j]);
                        //if the cost to the new city is cheaper, remember it
                        if (min > tempCost)
                        {
                            min = tempCost;
                            tempIndex = j;
                        }
                    }
                    //if the smallest route we found was infinity, we need to try to start in a different city
                    if (min == double.PositiveInfinity)
                    {
                        i = Cities.Length;
                        fail = true;
                    }
                    else //otherwise it is a valid route and we remember it
                    {
                        citiesVisited[tempIndex] = true;
                        route.Add(tempIndex);
                        currentCity = tempIndex;
                    }
                }
                //if we found a valid route, we can update the bssf otherwise we'll try a new city to start from
                if (!fail)
                    updateBssf(route);
            } while (costOfBssf() == double.PositiveInfinity || costOfBssf() == -1);


            //we then start the branch and a bound method
            //the first step is to make the matrix, this takes space complexity of O(n^2) and speed of O(n^2)
            double [,] matrix = new double[Cities.Length, Cities.Length];
            for(int i = 0; i < Cities.Length; i++)
            {
                for(int j = 0; j < Cities.Length; j++)
                {
                    if(i == j)
                    {
                        matrix[i,j] = double.PositiveInfinity;
                    }
                    else
                    {
                        matrix[i, j] = Cities[i].costToGetTo(Cities[j]);
                    }
                }
            }
            //we then reset our route, starting with route 0, and starting our citiesVisited matrix over
            route = new List<int>();
            route.Add(0);
            citiesVisited = new bool[Cities.Length];
            citiesVisited[0] = true;
            double lowerBound = 0;
            //we call the makeNewState method which take a matrix, reduces it and adds to the given lowerBound to make a state object
            WindowsFormsApplication1.State state = makeNewState(matrix, lowerBound, route, citiesVisited);
            //we also make our priority queue and add our new state we made
            WindowsFormsApplication1.PriorityQueue priorityQueue = new WindowsFormsApplication1.PriorityQueue();
            priorityQueue.add(state);

            //we initialize all of these variables to speed up our while loop
            int currentCityIndex;
            double tempLowerBound;
            double[,] tempMatrix;
            WindowsFormsApplication1.State tempState;
            List<int> tempRoute;
            bool[] tempCitiesVisited;
            double totalCost;
            int totalChildStates = 0;
            int totalPrunedStates = 0;
            int maxSizeQueue = 1;
            int tempQueueSize;

            //we will keep running through our while loop while there are states in the queue or until runs out
            //the longest method run in this loop is the for loop with O(n^3) and we could run for all states expanded which
            //could go up to O(n!) and space is how many states we expand on which could also go up to O(n!)
            while (priorityQueue.getSize() != 0)
            {
                //pop the state off the priority queue
                state = priorityQueue.deleteMin();
                //check if our lower bound is lower than our cost of bssf or we will skip it
                if (state.lowerBound >= costOfBssf())
                {
                    totalPrunedStates++;
                    continue;
                }
                //pull the variables from the state for convience
                route = state.route;
                citiesVisited = state.citiesVisited;
                lowerBound = state.lowerBound;
                matrix = state.reducedMatrix;
                //set our current city index to the last city in the route
                currentCityIndex = route[route.Count - 1];

                if(route.Count != Cities.Length)//route not finished
                {
                    //we expand on all cities were we haven't gone, this runs in O(n^3) time for the for loop and using the makeNewState method and O(n^2) space for the matrices
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        //don't expand if we've already gone to the city
                        if (citiesVisited[j]) 
                        {
                            continue;
                        }
                        totalChildStates++;
                        //add the cost to get to the new city and see if it goes over our upper bound, if so we can prune it
                        tempLowerBound = lowerBound + matrix[currentCityIndex, j];
                        if (tempLowerBound >= costOfBssf())
                        {
                            totalPrunedStates++;
                            continue;
                        }
                        //copy the reduced matrix from the previous city, infinity the row of the current city,
                        //infinity the column of the new city, deep copy the route adding the new city
                        //deep copy the cities visited, then from all that make a new state reducing the matrix down
                        //and updating the new lowerbound
                        tempMatrix = (double[,])matrix.Clone();
                        tempMatrix = infinityRow(tempMatrix, currentCityIndex, Cities.Length);
                        tempMatrix = infinityColumn(tempMatrix, j, Cities.Length);
                        tempRoute = route.ConvertAll(integer => integer);
                        tempRoute.Add(j);
                        tempCitiesVisited = (bool[])citiesVisited.Clone();
                        tempCitiesVisited[j] = true;
                        tempState = makeNewState(tempMatrix, tempLowerBound, tempRoute, tempCitiesVisited);
                        //if the new lower bound exceeds the bssf, prune it
                        if (tempState.lowerBound >= costOfBssf())
                        {
                            totalPrunedStates++;
                            continue;
                        }
                        //otherwise add the new state to the priority queue
                        priorityQueue.add(tempState);
                        tempQueueSize = priorityQueue.getSize();
                        if (maxSizeQueue < tempQueueSize)
                        {
                            maxSizeQueue = tempQueueSize;
                        }
                    }
                }
                else //finished the route
                {
                    //when the route is finished add on the value to travel back to city 0, and if we did better than
                    //the old bssf, update it
                    totalCost = lowerBound + matrix[currentCityIndex, 0];
                    if(totalCost < costOfBssf()) //found a better solution
                    {
                        updateBssf(route);
                        numSolutionsFound++;
                    }
                }

                //this checks if we have gone over our time limit and if we did then we'll break out
                if (timer.Elapsed > endTime)
                {
                    totalPrunedStates += priorityQueue.getSize();
                    break;
                }
                    
            }

            //update and return our results
            results[COST] = costOfBssf().ToString();                         
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = numSolutionsFound.ToString();

            return results;
        }

        private void updateBssf(List<int> route)
        {
            //takes the list of indexes and turns it into a list of cities
            Route = new ArrayList(Cities.Length);
            for(int i = 0; i < route.Count; i++)
            {
                Route.Add(Cities[route[i]]);

            }
            //updates the bssf
            bssf = new TSPSolution(Route);
            
        }

        private double[,] infinityRow(double[,] matrix, int index, int dimension)
        {
            //takes a given matrix and infinities the given index row
            for(int i = 0; i < dimension; i++)
            {
                matrix[index, i] = double.PositiveInfinity;
            }
            return matrix;
        }

        private double[,] infinityColumn(double[,] matrix, int index, int dimension)
        {
            //takes a given matrix and infinities the given index column
            for(int j = 0; j < dimension; j++)
            {
                matrix[j, index] = double.PositiveInfinity;
            }
            return matrix;
        }
        
       
        private WindowsFormsApplication1.State makeNewState(double[,] matrix, double lowerBound, List<int> route, bool[] citiesVisited)
        {
            //takes a matrix, reduces it updating and adding to the lower bound and returns a states
            //which contains the reduced matrix, lowerbound, route, and cities visited
            //the whole method runs in O(n^2) space and time
            double min;
            //runs through the all rows and finds the min and subtracts it off from all values in the row, adding the min to the lower bound
            //this runs in O(n^2) time for the two for loops and space for the matrix
            for (int i = 0; i < Cities.Length; i++)
            {
                min = double.PositiveInfinity;
                for(int j = 0; j < Cities.Length; j++)
                {
                    if(min > matrix[i, j])
                    {
                        min = matrix[i, j];
                    }
                }
                if (!double.IsPositiveInfinity(min))
                {
                    lowerBound += min;
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        matrix[i, j] -= min;
                    }
                }
                
            }
            //now it runs through all the columns doing the same, finding the min, subtracting it from all values in the column and adds it to the lowerbound
            //this runs in O(n^2) time and space
            for (int j = 0; j < Cities.Length; j++)
            {
                min = double.PositiveInfinity;
                for(int i = 0; i < Cities.Length; i++)
                {
                    if (min > matrix[i, j])
                    {
                        min = matrix[i, j];
                    }
                }
                if (!double.IsPositiveInfinity(min))
                {
                    lowerBound += min;
                    for(int i = 0; i < Cities.Length; i++)
                    {
                        matrix[i, j] -= min;
                    }
                }
            }

            //after reducing the matrix and updating the lowerBound, it makes a new state and returns it
            return new WindowsFormsApplication1.State(matrix, lowerBound, route, citiesVisited);
        }
        

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            string[] results = new string[3];

            Stopwatch greedyTimer = new Stopwatch();
            greedyTimer.Start();
            bssf = GreedyImplementation(0);
            greedyTimer.Stop();

            results[COST] = "" + bssf.costOfRoute();    // load results into array here, replacing these dummy values
            results[TIME] = greedyTimer.Elapsed.ToString();
            results[COUNT] = "1";

            return results;
        }

        private TSPSolution GreedyImplementation(int start)
        {
            //find our bssf and start at node 0
            //Uses a greedy algorithm that takes the closest node and continues on.
            //If in the hard version we end up getting stuck with an infinite route, we start on a different node
            //and run the greedy algorithm again.
            //If none of the nodes produces a viable solution, we run the default algoritm which finds a random route.
            List<int> route;
            //min keeps track of the smallest route we've found so far
            double min;
            //already gone keeps track of what cities we have already gone to to not look at those again
            bool[] citiesVisited;
            double tempCost;
            int tempIndex;
            int currentCity = -1;
            //if fail is true then we found an infinite route and nee
            bool fail;
            do //it is unlikely but we could have to run through all cities, and then go to the default solution so this could take speed O(n^3) with space complexity of O(n) for arrays and list
            {

                fail = false;
                currentCity++;
                citiesVisited = new bool[Cities.Length];
                citiesVisited[currentCity] = true;
                route = new List<int>();
                route.Add(currentCity);

                if (currentCity >= Cities.Length) //if all greedy attempts don't solve it, use the default solution
                {
                    defaultSolveProblem();
                    break;
                }
                //begin at city start and find a route to all cities
                //we connect to every city and run through this n times which contains a O(n) method so it is O(n^2) with the space complexity of O(n) for arrays and list
                for (int i = start; i < Cities.Length - 1; i++)
                {
                    min = double.PositiveInfinity;
                    tempIndex = -1;
                    //look at the route from our currentCity all available cities
                    //on average, we will run through this for loop half the time or O(n/2) which is O(n) speed complexity with the same size complexity (for cities visited)
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        //skip the routes where we have already gone
                        if (citiesVisited[j])
                            continue;
                        tempCost = Cities[currentCity].costToGetTo(Cities[j]);
                        //if the cost to the new city is cheaper, remember it
                        if (min > tempCost)
                        {
                            min = tempCost;
                            tempIndex = j;
                        }
                    }
                    //if the smallest route we found was infinity, we need to try to start in a different city
                    if (min == double.PositiveInfinity)
                    {
                        i = Cities.Length;
                        fail = true;
                    }
                    else //otherwise it is a valid route and we remember it
                    {
                        citiesVisited[tempIndex] = true;
                        route.Add(tempIndex);
                        currentCity = tempIndex;
                    }
                }


                for (int i = 0; i < start; i++)
                {
                    min = double.PositiveInfinity;
                    tempIndex = -1;
                    //look at the route from our currentCity all available cities
                    //on average, we will run through this for loop half the time or O(n/2) which is O(n) speed complexity with the same size complexity (for cities visited)
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        //skip the routes where we have already gone
                        if (citiesVisited[j])
                            continue;
                        tempCost = Cities[currentCity].costToGetTo(Cities[j]);
                        //if the cost to the new city is cheaper, remember it
                        if (min > tempCost)
                        {
                            min = tempCost;
                            tempIndex = j;
                        }
                    }
                    //if the smallest route we found was infinity, we need to try to start in a different city
                    if (min == double.PositiveInfinity)
                    {
                        i = Cities.Length;
                        fail = true;
                    }
                    else //otherwise it is a valid route and we remember it
                    {
                        citiesVisited[tempIndex] = true;
                        route.Add(tempIndex);
                        currentCity = tempIndex;
                    }
                }
                //if we found a valid route, we can update the bssf otherwise we'll try a new city to start from
                if (!fail)
                    updateBssf(route);
            } while (costOfBssf() == double.PositiveInfinity || costOfBssf() == -1);
            return bssf;
        }

        public string[] fancySolveProblem()
        {
            string[] results = new string[3];

            int population = 10000;
            int randomnums = 10;
            List<Link> GlobalBest = new List<Link>();
            double bestDistance = double.PositiveInfinity;

            List<List<Link>> routes = new List<List<Link>>();

            for(int i = 0; i < Cities.Length; i++)
            {
                TSPSolution tsp = GreedyImplementation(i);
                if((tsp.costOfRoute()) < bestDistance)
                {
                    GlobalBest = Linkify(tsp);
                    bestDistance = tsp.costOfRoute();
                }
                routes.Add(Linkify(GreedyImplementation(i)));
            }

            //Making a bad assumption that default won't return the best solution because that requires rewriting the defaul method
            for(int i = routes.Count; i< population; i++)
            {
                defaultSolveProblem();
                routes.Add(Linkify(bssf));
            }

            Random rand = new Random();
            Stopwatch fancyTimer = new Stopwatch();
            fancyTimer.Start();
            while(fancyTimer.Elapsed.Milliseconds < (time_limit*1000))
            {
                List<int> GrabRandom = new List<int>();
                List<double> distances = new List<double>();

                for (int i = 0; i < randomnums; i++)
                {
                    GrabRandom.Add(rand.Next(population));
                    List<Link> ls = routes[GrabRandom[i]];
                    double distance = 0;
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        distance += Cities[ls[j].first].costToGetTo(Cities[ls[j].second]);
                    }
                    distances.Add(distance);
                }

                //sorting here because copying the GrabRandom array and giant routes array was ridiculous


                List<Link> child = crossover(routes[GrabRandom[0]], routes[GrabRandom[1]]);
                List<Link> worst = routes[GrabRandom[GrabRandom.Count - 1]];
                double childDistance = 0.0;
                double worstDistance = 0.0;
                for(int i = 0; i < Cities.Length; i++)
                {
                    childDistance += Cities[child[i].first].costToGetTo(Cities[child[i].second]);
                    worstDistance += Cities[worst[i].first].costToGetTo(Cities[worst[i].second]);
                }

                if(childDistance < worstDistance)
                {
                    routes[GrabRandom[GrabRandom.Count - 1]] = child;
                    if (childDistance < bestDistance)
                    {
                        GlobalBest = child;
                        bestDistance = childDistance;
                    }
                }

                
            }

            results[COST] = "" + bestDistance;    // load results into array here, replacing these dummy values
            results[TIME] = "60.0";
            results[COUNT] = "1";

            return results;
        }

        //Take a tspsolution and make it a list of links
        private List<Link> Linkify(TSPSolution tsp)
        {
            List<Link> route = new List<Link>();
            City first, second; //Route is an Arraylist of objects and setting them as cities in a single line as a new link was a little crazy

            //Loop through every city and make links from the two cities that are connected
            for(int i =0; i < tsp.Route.Count - 1; i++)
            {
                first = (City)tsp.Route[i];
                second = (City)tsp.Route[i + 1];
                route.Add(new Link(first.position, second.position));
            }

            first = (City)tsp.Route[tsp.Route.Count - 1];
            second = (City)tsp.Route[0];
            route.Add(new Link(first.position, second.position)); //Don't forget the loop around

            return route;
        }
        
        
        #endregion
    }

}
