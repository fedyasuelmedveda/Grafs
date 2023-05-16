using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using OpenTK;
using OpenTK.Graphics.OpenGL;
using DifferentialEquations;
using System.Windows.Forms.VisualStyles;

namespace Graphs
{
    public partial class Form1 : Form
    {
        Graph g;
        int grid = 100;
        int vertexForComponent = 0;
        public Form1()
        {
            InitializeComponent();
        }
        void Vertexes()
        {
            List<DoublePoint> v = g.GetVertexes();
            int grid = 20;
            for(int j = 0; j < v.Count; j++)
            {
                if (mark[j] == true)
                {
                    GL.Begin(PrimitiveType.TriangleFan);
                    GL.Color3(1.0f, 0.0f, 0.0f);
                    double delta = 2.0 / grid;
                    DoublePoint vertex = v[j];
                    for (int i = 0; i < grid; i++)
                    {
                        GL.Vertex2(vertex.X + 0.2 * Math.Cos(Math.PI * delta * i), vertex.Y + 0.2 * Math.Sin(Math.PI * delta * i));
                    }
                    GL.End();
                    Number(j, vertex);
                }
            }
        }

        void Edges()
        {
            
            List<List<double>> e = g.GetEdges();

            foreach(List<double> edge in e)
            {
                Edge(edge);
            }
        }
        void Number(int n,DoublePoint location)
        {
            GL.Begin(PrimitiveType.Lines);
            GL.Color3(0.0f, 0.0f, 1.0f);
            double size = 0.05;
            if (n == 0)
            {
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);

                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X -1 * size, location.Y - 2 * size);
            }

            if (n == 1)
            {
                GL.Vertex2(location.X , location.Y - 2 * size);
                GL.Vertex2(location.X, location.Y + 2 * size);
            }
            if (n == 2)
            {
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y);

                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y);
            }
            if (n == 3)
            {
                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y);
            }
            if (n == 8)
            {
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X + 1 * size, location.Y - 2 *size);

                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

            }
            if (n == 7)
            {
                //GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);
                //GL.Vertex2(location.X - 1 * size, location.Y);

                //GL.Vertex2(location.X - 1 * size, location.Y);
                //GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);

                //GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                //GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);

                //GL.Vertex2(location.X + 1 * size, location.Y);
                //GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

            }

            if (n == 6)
            {
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);

                //GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                //GL.Vertex2(location.X + 1 * size, location.Y);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

            }

            if (n == 5)
            {
                //GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);
                //GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);

                //GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                //GL.Vertex2(location.X + 1 * size, location.Y);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

            }
            if (n == 4)
            {
                //GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);
                //GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);

                //GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                //GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y);

               //GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                //GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

            }
            if (n == 9)
            {
                //GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);
                //GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y - 2 * size);
                GL.Vertex2(location.X - 1 * size, location.Y - 2 * size);

                GL.Vertex2(location.X + 1 * size, location.Y);
                GL.Vertex2(location.X - 1 * size, location.Y);

                GL.Vertex2(location.X - 1 * size, location.Y + 2 * size);
                GL.Vertex2(location.X + 1 * size, location.Y + 2 * size);

            }

            GL.End();
        }
        void DrawGraph()
        {
            GL.Clear(ClearBufferMask.ColorBufferBit);
            Edges();
            Vertexes();
            
            glControl1.SwapBuffers();
        }
        int Cnk(int k, int n)
        {
            int res = 1;
            for (int i = n; i > n - k; i--)
            {
                res *= i;
            }
            for (int i = 1; i <= k; i++)
            {
                res /= i;
            }
            return res;
        }
        void Edge(List<double> edges)
        {
            
            List<DoublePoint> vertexes = g.GetVertexes();
            Vector3[,] colors = g.GetColorMatrix();
            DoublePoint u = vertexes[(int)edges[0]];
            DoublePoint v = vertexes[(int)edges[1]];
            DoublePoint m = (1.0/ 2) * (u + v);
            DoublePoint ortv = new DoublePoint(u.Y - v.Y, v.X - u.X);
            ortv = (0.5 / ortv.Norm()) * ortv;
            //ortv = 0.1 * ortv;
            if (mark[(int)edges[0]] == true && mark[(int)edges[1]] == true){
                for (int k = 2; k < edges.Count; k++)
                {
                    double arrowx = 0;
                    double arrowy = 0;
                    GL.Begin(PrimitiveType.LineStrip);
                    GL.Color3(colors[(int)edges[0], (int)edges[1]]);
                    DoublePoint[] vec = new DoublePoint[] { u, m + (k - 1) * ortv, v };
                    double delta = 1 / (double)grid;
                    for (int i = 0; i <= grid; i += 1)
                    {
                        double x = 0;
                        double y = 0;
                        for (int j = 0; j < 3; j++)
                        {
                            double product = Cnk(j, 2) * Math.Pow(1 - i * delta, 2 - j) * Math.Pow(i * delta, j);
                            x += (product * vec[j].X);
                            y += (product * vec[j].Y);

                        }
                        GL.Vertex2(x, y);
                        if (i == grid / 2)
                        {
                            arrowx = x;
                            arrowy = y;
                        }

                    }
                    GL.End();
                    GL.Begin(PrimitiveType.Triangles);
                    GL.Color3(0.0f, 1.0f, 0.0f);
                    double n = 0.1 / (v - u).Norm();
                    GL.Vertex2(arrowx, arrowy);
                    GL.Vertex2(arrowx - n * (v.X - u.X) + 0.1 * ortv.X, arrowy - n * (v.Y - u.Y) + 0.1 * ortv.Y);
                    GL.Vertex2(arrowx - n * (v.X - u.X) - 0.1 * ortv.X, arrowy - n * (v.Y - u.Y) - 0.1 * ortv.Y);
                    GL.End();

                }
            }
        }

        void Shortest()
        {
            List<List<double>> e = g.GetEdges();
            List<List<double>> newedges = new List<List<double>>();
            double[,] newmatrix = g.GetAdjacencyMatrix();
            foreach (List<double> edge in e)
            {
                double min = 10000;
                for(int i = 2; i < edge.Count; i++)
                {
                    if (edge[i]<min)
                        min = edge[i];
                }
                newedges.Add( new List<double> { edge[0], edge[1], min });
                newmatrix[(int)edge[0], (int)edge[1]] = min;
                richTextBox1.Text += $"{edge[0]} -> {edge[1]}:";
                richTextBox1.Text += '\n';

                richTextBox1.Text += min.ToString();
                richTextBox1.Text += '\n';
            }
            richTextBox1.Text += '\n';
            g.SetEdges(newedges);
            g.SetAdjacencyMatrix(newmatrix);

        }
        private void glControl1_Load(object sender, EventArgs e)
        {
            GL.ClearColor(1.0f, 1.0f, 1.0f, 1.0f);
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            GL.Ortho(-2, 2, -2, 2, -2, 2);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();
            double[,] matrix = new double[,] { {0,1,0,1 },{1,0,1,0 },{0,1,0,1 },{ 1,0,1,0} };
            for(int i =0;i < 20; i++)
            {
                mark[i] = true;
            }
            
        }

        private void button1_Click(object sender, EventArgs e)
        {
            g = new Graph('r', 11);
            
            for (int i = 0; i < 20; i++)
            {
                mark[i] = true;
            }
            DrawGraph();
        }

        private void button2_Click(object sender, EventArgs e)
        {
            Shortest();
            DrawGraph();
        }
        DoublePoint GetCoords(DoublePoint p)
        {
            return new DoublePoint((double)p.X / glControl1.Size.Width * 4 - 2, 2 - (double)p.Y / glControl1.Size.Height * 4);
        }
        public DoublePoint Convert(Point p)
        {
            return new DoublePoint((double)p.X / glControl1.Size.Width * 4 - 2, 2 - (double)p.Y / glControl1.Size.Height * 4);

        }

        bool[] mark = new bool[20]; 
        void DfsComponent(int v, List<List<double>> neighbours)
        {
            mark[v] = true;
            for(int i = 0; i < neighbours[v].Count; i++) 
            {

                    int w = (int)neighbours[v][i];
                    if (!mark[w])
                        DfsComponent(w, neighbours);
            }
        }
        private void Component()
        {
            for (int i = 0; i < 20; i++)
            {
                mark[i] = false;
            }
            List<List<double>> neighbours = g.GetNeighbours();
            DfsComponent(vertexForComponent, neighbours);
        }
        private void glControl1_MouseClick(object sender, MouseEventArgs e)
        {
            List<DoublePoint> vertexes = g.GetVertexes();
            DoublePoint m = Convert(e.Location);
            for (int i = 0; i < vertexes.Count; i++)
            {
                double x = vertexes[i].X;
                double y = vertexes[i].Y;
                
                if (Math.Abs(x - m.X) < 0.1 && Math.Abs(y - m.Y) < 0.1)
                {
                    vertexForComponent = i;
                    //richTextBox1.Text = m.X.ToString() + '\n' + m.Y.ToString() + '\n' + x.ToString() + '\n' + y.ToString();

                }
            }
            //Component();
            double[,] ways = g.BellmanFord(vertexForComponent);
            int size = vertexes.Count;
            //richTextBox1.Clear();
            richTextBox1.Text += "WeightsOfWays:";
            richTextBox1.Text += '\n';
            for (int i = 0; i < size; i++)
            {
                richTextBox1.Text += ways[size - 1, i].ToString();
                richTextBox1.Text += ' ';
            }
            richTextBox1.Text += '\n';
            richTextBox1.Text += "WayTo1:";
            richTextBox1.Text += '\n';
            double[,] vertexesInWay = g.GetI();
            
            int currentVertex = 1;
            int newVertex = 1;
            for (int i = size - 1; i >= 0; i--)
            {
                if (newVertex != currentVertex)
                {
                    richTextBox1.Text += currentVertex;
                    richTextBox1.Text += ' ';
                }
                currentVertex = newVertex;
                if (currentVertex != -1)
                    newVertex = (int)vertexesInWay[i, currentVertex];

            }

            richTextBox1.Text += vertexForComponent;
            richTextBox1.Text += '\n';
            richTextBox1.Text += "WeightsOfWays:";
            richTextBox1.Text += '\n';
            double[] waysShort = g.BellmanFordShort(vertexForComponent);
            //richTextBox1.Clear();
            for (int i = 0; i < size; i++)
            {
                richTextBox1.Text += waysShort[i].ToString();
                richTextBox1.Text += ' ';
            }
            richTextBox1.Text += '\n';
            richTextBox1.Text += "WayTo1";
            richTextBox1.Text += '\n';
            double[] vertexesInWayShort = g.GetIShort();
            //richTextBox1.Text += '\n';
            currentVertex = 1;
            newVertex = 1;
            for (int i = size - 1; i >= 0; i--)
            {
                if (newVertex != currentVertex)
                {
                    richTextBox1.Text += currentVertex;
                    richTextBox1.Text += ' ';
                }
                currentVertex = newVertex;
                if (currentVertex != -1)
                {
                    if ((int)vertexesInWayShort[currentVertex] != -1)
                    {
                        g.SetColorValue((int)vertexesInWayShort[currentVertex], currentVertex, new Vector3(1.0f, 0.0f, 0.0f));
                    }
                    newVertex = (int)vertexesInWayShort[currentVertex];
                }
            }
            richTextBox1.Text += vertexForComponent;
            richTextBox1.Text += '\n';

            richTextBox1.Text += "WeightsOfWays:";
            richTextBox1.Text += '\n';
            waysShort = g.Dejkstra(vertexForComponent);
            //richTextBox1.Clear();
            for (int i = 0; i < size; i++)
            {
                richTextBox1.Text += waysShort[i].ToString();
                richTextBox1.Text += ' ';
            }
            richTextBox1.Text += '\n';
            richTextBox1.Text += "WayTo1";
            richTextBox1.Text += '\n';
            vertexesInWayShort = g.GetIShort();
            currentVertex = 1;
            newVertex = 1;
            for (int i = size - 1; i >= 0; i--)
            {
                if (newVertex != currentVertex)
                {
                    richTextBox1.Text += currentVertex;
                    richTextBox1.Text += ' ';
                }
                currentVertex = newVertex;
                if (currentVertex != -1)
                {
                    if ((int)vertexesInWayShort[currentVertex] != -1)
                    {
                        g.SetColorValue((int)vertexesInWayShort[currentVertex], currentVertex, new Vector3(1.0f, 0.0f, 0.0f));
                    }
                    newVertex = (int)vertexesInWayShort[currentVertex];
                }
            }
            richTextBox1.Text += vertexForComponent;
            richTextBox1.Text += '\n';
            g.ResetWays();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            DrawGraph();
        }
    }
}
