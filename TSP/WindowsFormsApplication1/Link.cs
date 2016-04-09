using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP
{
    class Link
    {
        private int _first, _second;
        public Link(int x, int y)
        {
            _first = x;
            _second = y;
        }

        public int first
        {
            get{return _first; }
            set {_first = value; }
        }

        public int second
        {
            get { return _second; }
            set { _second = value; }
        }

        //returns a 0 if they are the same.
        public int compareTo(Link next)
        {
            if (next.first == _first)
                if (next.second == _second)
                    return 0;
            return 1;
        }

        public string toString()
        {
            return first + ", " + second;
        }
    }
}
