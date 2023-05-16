using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using DifferentialEquations;
using OpenTK;
using OpenTK.Graphics.ES11;

namespace Graphs
{
    internal class Graph
    {
        private List<DoublePoint> Vertexes;
        private List<List<double>> Edges;
        private List<List<double>> Neighbours;
        private double[,] ways;
        private double[,] weightMatrix;
        private double[,] adjacencyMatrix;
        private Vector3[,] colorMatrix;
        double[,] I;
        double[] IShort;
        double[] waysShort;
        Random random = new Random();
        int size;
        bool mark;

        public Graph(double[,] matrix)
        {
            
            adjacencyMatrix = matrix;
            size = matrix.GetLength(0);
            ways = new double[size, size];
            colorMatrix = new Vector3[size, size];
            ResetWays();
            Vertexes = new List<DoublePoint>();
            VertexesCoordinates();
            Edges = new List<List<double>>();
            Neighbours = new List<List<double>>(20);
            for(int i = 0; i < size; i++)
            {
                for(int j= 0; j < size; j++)
                {
                    if (i != j)
                    {
                        if (adjacencyMatrix[i, j] > 0)
                        {
                            Edges.Add(new List<double> { i, j });
                            for (int k = 0; k < adjacencyMatrix[i, j]; k++)
                            {
                                int weight = random.Next(20);
                                Edges[Edges.Count - 1].Add(weight);
                            }
                            Neighbours[(int)i].Add(j);
                            colorMatrix[i,j] = new Vector3(0.0f, 1.0f, 0.0f);
                        }
                    }
                }
            }
        }
        public Graph(char c,int s)
        {
            if (c == 'e')
            {
                Vertexes = new List<DoublePoint>();
                Edges = new List<List<double>>();
                size = s;
                for(int i=0; i < size; i++)
                    for (int j = 0; j < size; j++)
                        adjacencyMatrix[i, j] = 0;
            }
            if (c == 'r')
            {
                size = s;
                ways = new double[size, size];
                ResetWays();
                random = new Random();
                Vertexes = new List<DoublePoint>();
                Edges = new List<List<double>>();
                colorMatrix = new Vector3[size, size];
                Neighbours = new List<List<double>>();
                for(int i = 0; i< 20; i++)
                {
                    Neighbours.Add(new List<double>());
                }
                
                adjacencyMatrix = new double[size, size];
                VertexesCoordinates();
                for (int i = 0; i < size; i++)
                    for (int j = 0; j < size; j++)
                    {
                        if (i != j) { 
                        int r = random.Next(2);
                        adjacencyMatrix[i, j] = r;
                            //int r = 1;
                            if (r > 0)
                            {
                                Edges.Add(new List<double> { i, j });
                                for (int k = 0; k < adjacencyMatrix[i, j]; k++)
                                {
                                    int weight = random.Next(1,20);
                                    Edges[Edges.Count - 1].Add(weight);
                                }
                                Neighbours[(int)i].Add(j); 
                                colorMatrix[i, j] = new Vector3(0.0f, 1.0f, 0.0f);

                            }
                            else if(r==0)
                                adjacencyMatrix[(int)i, j] = double.PositiveInfinity;
                        }
                    }
            }/*
            if (c == 's')
            {
                Vertexes = new List<DoublePoint>();
                Edges = new List<List<double>>();
                size = s;
                VertexesCoordinates();
                for (int i = 0; i < size; i++)
                    for (int j = i + 1; j < size; j++)
                    {
                        int r = random.Next(1);
                        adjacencyMatrix[i, j] = r;
                        if (r == 1)
                        {
                            Edges.Add(new DoublePoint(i, j));
                            Edges.Add(new DoublePoint(j, i));
                        }
                    }
            }*/
        }
        public List<DoublePoint> GetVertexes()
        {
            return Vertexes;
        }
        public List<List<double>> GetEdges()
        {
            return Edges;
        }
        public void SetEdges(List<List<double>> newEdges)
        {
            Edges = newEdges;
        }
        public double [,] GetAdjacencyMatrix()
        {
            return adjacencyMatrix;
        }
        public Vector3[,] GetColorMatrix()
        {
            return colorMatrix;
        }
        public void SetColorValue(int i, int j, Vector3 v) {
            colorMatrix[i, j] = v;
        }
        public void SetAdjacencyMatrix(double[,] newMatrix)
        {
            adjacencyMatrix = newMatrix;
        }
        public List<List<double>> GetNeighbours()
        {
            return Neighbours;
        }
        private void VertexesCoordinates()
        {
            double delta = Math.PI * 2 / size;
            for(int i = 0; i < size; i++)
            {
                Vertexes.Add(new DoublePoint(1.5*Math.Cos(i*delta), 1.5*Math.Sin(i * delta)));
            }
        }
        public void ResetWays()
        {
            for(int i = 0; i < size; i++)
            {
                for(int j = 0; j < size; j++)
                {
                    ways[i, j] = double.PositiveInfinity;
                }
            }
        }
        public double[,] GetI()
        {
            return I;
        }
        public double[] GetIShort()
        {
            return IShort;
        }

        public double[,] BellmanFord(int s)
        {
            I = new double[size, size];
            //IShort = new double[size];
            
            for (int i = 0; i < size; i++)
            {
                ways[0,i] = double.PositiveInfinity;
                ways[1, i] = adjacencyMatrix[s, i];
                if (ways[1, i] < double.PositiveInfinity)
                {
                    I[1, i] = s;
                    //IShort[i] = s;
                }
                else
                {
                    I[1, i] = -1;
                    //IShort[i] = -1;
                }
            }
            ways[0, s] = 0;
            for (int i = 2; i < size; i++)
            {
                bool mark = false;
                for (int j = 0; j < size; j++)
                {
                    double min = ways[i - 1, j];
                    int minvertex = j;
                    for (int k = 0; k < size; k++)
                    {
                        if(ways[i - 1, k] + adjacencyMatrix[k, j]< min)
                        {
                            min = ways[i - 1, k] + adjacencyMatrix[k, j];
                            minvertex = k;
                            mark = true;

                        }

                    }
                    ways[i, j] = min;
                    I[i, j] = minvertex;
                    //IShort[j] = minvertex;
                }
                if (!mark)
                {
                    for (int j = i+1; j < size; j++)
                    {
                        for(int k = 0; k < size; k++)
                        {
                            ways[j, k] = ways[i, k];
                            I[j, k] = k;
                        }
                        
                    }
                    break;
                }
            }
            return ways;
        }

        public double[] BellmanFordShort(int s)
        {
            //I = new double[size, size];
            IShort = new double[size];
            waysShort = new double[size];
            waysShort[s] = s;
            for (int i = 0; i < size; i++)
            {
                waysShort[i] = adjacencyMatrix[s, i];
                if (waysShort[i] < double.PositiveInfinity)
                {
                    IShort[i] = s;
                }
                else
                {
                    IShort[i] = -1;
                }
            }
            for (int i = 2; i < size; i++)
            {
                bool mark1 = true;
                for (int j = 0; j < size; j++)
                {
                    double min = waysShort[j];
                    int minvertex = j;
                    bool mark2 = false;
                    for (int k = 0; k < size; k++)
                    {
                        if (waysShort[k] + adjacencyMatrix[k, j] < min)
                        {
                            min = waysShort[k] + adjacencyMatrix[k, j];
                            minvertex = k;
                            mark2 = true;
                            mark1 = false;
                        }

                    }
                    waysShort[j] = min;
                    if(mark2)
                    IShort[j] = minvertex;
                }
                if (mark1)
                    break;
            }
            return waysShort;
        }
        public double[] Dejkstra(int s) {
            waysShort = new double[size];
            IShort = new double[size];
            bool[] mark = new bool[size];
            for(int i = 0;i < size; i++)
            {
                mark[i] = false;
                waysShort[i] = double.PositiveInfinity;
            }
            waysShort[s] = 0;
            int current = s;
            while (true)
            {
                mark[current] = true;
                for(int i =0; i< Neighbours[current].Count; i++)
                {
                    
                    int t = (int)Neighbours[current][i];
                    if (!mark[t])
                        if (waysShort[t] > waysShort[current] + adjacencyMatrix[current, t])
                        {
                            waysShort[t] = waysShort[current] + adjacencyMatrix[current, t];
                            IShort[t] = current;
                        }
                }
                bool ifStop = false;
                bool ifStop1 = false;
                double min = double.PositiveInfinity;
                int minVertex = current;
                for(int i = 0; i < size; i++)
                {
                    if (!mark[i])
                    {
                        ifStop1 = true;
                        if (waysShort[i] < double.PositiveInfinity)
                        {
                            ifStop = true;
                            if (waysShort[i] < min)
                            {
                                min = waysShort[i];
                                minVertex = i;
                            }

                        }
                    }
                }
                if (!ifStop || !ifStop1)
                    break;
                current = minVertex;

            }
            return waysShort;
        }
    }
}
