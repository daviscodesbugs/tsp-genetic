using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WindowsFormsApplication1
{
    class PriorityQueue
    {
        //the heap queue has the heap array
        //it orders states in the heap based on their returned priority value
        //this heap queue has a space complexity of O(n) for the heap arrays
        private State[] heap;
        private int size;
       

        //returns the current number of valid values in the queue
        //this runs in O(1) time
        public int getSize()
        {
            return size;
        }

        //make queue makes the queue initializing the potential size of the queue to 100 which can grow as needed
        public PriorityQueue()
        {
            size = 0;
            heap = new State[100];
        }

        //adds a state to the queue to the bottom
        //before it adds the state it makes sure there is room, if not, it doubles its size, adds the state to the bottom
        //of the heap and then checks if the parent is bigger in which case it will swap the two and continue looking at the
        //parent
        //add to the heap, it takes O(log(n)) to reorganize the heap, and if it needs to be resized it can add to this
        public void add(State newState)
        {
            if(size == heap.Length)
            {
                Array.Resize<State>(ref heap, size * 2);
            }
            heap[size] = newState;
            size++;
            checkParentSwap(size - 1);
            
        }
        

        //the delete min method pulls out the smallest node and replaces it with the last node
        //and then runs a check child swap method which recursively switches the node with the smaller of 
        //its children.
        //this method also adjusts the heap position array to keep track of the indeces in the heap
        //this method runs in O(log(n)) time as it may call checkChildSwap that many times
        public State deleteMin()
        {
            State minState = heap[0];
            heap[0] = heap[size - 1];
            heap[size - 1] = null;
            size--;
            if (heap[0] != null)
                checkChildSwap(0);
            return minState;
        }

        //in the heap array (i-1)/2 represents the parent in the heap and then checks to see if there is a parent
        //and if there is it checks if the parent is bigger, if so it will swap the two and continue to check if the
        //new parent is bigger
        //this runs in O(1) times since it is just comparing, but as it runs recursively it could call itself O(log(n)) times
        private void checkParentSwap(int index)
        {
            int parentIndex = (index - 1) / 2;
            State parent = heap[parentIndex];
            if (parent != null) //there is a parent
            {
                if (parent.getPriorityValue() > heap[index].getPriorityValue())
                {
                    swap(index, parentIndex);
                    checkParentSwap(index);
                }
            }
        }

        //this method does a similar thing to the checkParentSwap method except that it compares the two children
        //and sees which is smaller and also has to verify that both children may or may not exist
        //so it has a higher complexity since there are so many possibilities
        //this runs in O(1) time but it may call itself O(log(n)) times as it recurses
        private void checkChildSwap(int index)
        {
            int child1Index = 2 * index + 1;
            int child2Index = 2 * index + 2;
            
            
            if (child1Index < size) //child 1 exists
            {
                State child1 = heap[child1Index];
                if (child2Index < size) //child 2 exists
                {
                    State child2 = heap[child2Index];
                    if (child1.getPriorityValue() < child2.getPriorityValue()) //child1 is smaller than child2
                    {
                        if (child1.getPriorityValue() < heap[index].getPriorityValue()) //swap child1 if smaller
                        {
                            swap(index, child1Index);
                            checkChildSwap(index);
                        }
                    }
                    else //child2 is smaller than child1
                    {
                        if (child2.getPriorityValue() < heap[index].getPriorityValue()) //swap child2 if smaller
                        {
                            swap(index, child2Index);
                            checkChildSwap(index);
                        }
                    }
                }
                else // only 1 child exists
                {
                    if (child1.getPriorityValue() < heap[index].getPriorityValue()) //swap child 1 if smaller
                    {
                        swap(index, child1Index);
                        checkChildSwap(index);
                    }
                }
            }
        }

        //this method swaps the value of two states in the heap
        //arrays
        //this runs in O(1) time to switch the nodes' positions
        private void swap(int index1, int index2)
        {
            State tempState = heap[index1];
            heap[index1] = heap[index2];
            heap[index2] = tempState;

        }
    
    }
}
