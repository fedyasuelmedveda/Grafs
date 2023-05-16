using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DifferentialEquations
{
    public struct DoublePoint
    {
        public double X;
        public double Y;
        public DoublePoint(double x1, double x2)
        {
            X = x1;
            Y = x2;
        }
        public DoublePoint(Point p)
        {
            X = p.X;
            Y = p.Y;
        }
        public double Norm()
        {
            return Math.Sqrt(X*X + Y*Y);
        }
        public static DoublePoint operator +(DoublePoint x, DoublePoint y)
        {
            return new DoublePoint(x.X + y.X, x.Y + y.Y);
        }
        public static DoublePoint operator -(DoublePoint x, DoublePoint y)
        {
            return new DoublePoint(x.X - y.X, x.Y - y.Y);
        }

        public static DoublePoint operator *(double h, DoublePoint x)
        {
            return new DoublePoint(h * x.X, h * x.Y);
        }
    }

}
