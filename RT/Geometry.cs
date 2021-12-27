using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.IO;
using System.Drawing;
using System.Globalization;

using static RT.Athens;

namespace RT
{
    static class Geometry
    {
        public static Scene scene;

        public static double Distance(Point3D p1, Point3D p2)
        {
            return Math.Sqrt(Math.Pow(p1.X - p2.X, 2) + Math.Pow(p1.Y - p2.Y, 2) + Math.Pow(p1.Z - p2.Z, 2));
        }

        public static double dotProduct(Point3D p1, Point3D p2)
        {
            return p1.X * p2.X + p1.Y * p2.Y + p1.Z * p2.Z;
        }

        public static Point3D normalize(Point3D p)
        {
            double d = Distance(new Point3D(0, 0, 0), p);
            return new Point3D(p.X / d, p.Y / d, p.Z / d);
        }

        //Отражение
        public static Point3D reflect(Point3D I, Point3D N)
        {
            double dot1 = dotProduct(I, N);
            return new Point3D(I.X - N.X * 2.0 * dot1, I.Y - N.Y * 2.0 * dot1, I.Z - N.Z * 2.0 * dot1);
        }

        //Проломление
        public static Point3D refract(Point3D dir, Point3D norm, double ref_index)
        {
            double dot1 = -1 * (dir * norm);
            double oldInd = 1;
            double new_index = ref_index;
            Point3D n = norm;
            if (dot1 < 0)
            {
                dot1 = -dot1;
                double t = oldInd;
                oldInd = new_index;
                new_index = t;
                n = norm * (-1);
            }
            double eta = oldInd / new_index;
            double k = 1 - eta * eta * (1 - dot1 * dot1);
            return k < 0 ? new Point3D(0, 0, 0) : dir * eta + n * (eta * dot1 - Math.Sqrt(k));
        }

        public static double PointPlaneDistance(Polygon pol, Point3D p, bool mod = true)
        {
            Point3D p1 = pol.points[0];
            Point3D p2 = pol.points[1];
            Point3D p3 = pol.points[2];

            double x1 = p1.X, y1 = p1.Y, z1 = p1.Z;
            double x2 = p2.X, y2 = p2.Y, z2 = p2.Z;
            double x3 = p3.X, y3 = p3.Y, z3 = p3.Z;

            double a = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
            double b = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
            double c = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
            double d = (-1) * (a * x1 + b * y1 + c * z1);
            if (mod)
            {
                return Math.Abs(a * p.X + b * p.Y + c * p.Z + d) / Math.Sqrt(a * a + b * b + c * c);
            }
            else
            {
                return (a * p.X + b * p.Y + c * p.Z + d) / Math.Sqrt(a * a + b * b + c * c);
            } 
        }

        public static Point3D CreateNormal(Polygon polygon, bool clockwise)
        {
            Point3D v1 = new Point3D(polygon.points[0].X - polygon.points[1].X,
                                     polygon.points[0].Y - polygon.points[1].Y,
                                     polygon.points[0].Z - polygon.points[1].Z);
            Point3D v2 = new Point3D(polygon.points[2].X - polygon.points[1].X,
                                     polygon.points[2].Y - polygon.points[1].Y,
                                     polygon.points[2].Z - polygon.points[1].Z);
            Point3D normalv = new Point3D();

            if (clockwise)
            {
                normalv = new Point3D(v1.Z * v2.Y - v1.Y * v2.Z, v1.X * v2.Z - v1.Z * v2.X, v1.Y * v2.X - v1.X * v2.Y);
            }
            else
            {
                normalv = new Point3D(v1.Y * v2.Z - v1.Z * v2.Y, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);
            }

            double dist = Distance(new Point3D(0, 0, 0), normalv);
            return new Point3D(normalv.X / dist, normalv.Y / dist, normalv.Z / dist);
        }

        public static double CalcAngleSum(Point3D q, Polygon p)
        {
            int n = p.points.Count;
            double m1, m2;
            double anglesum = 0, costheta = 0;
            double epsilon = 0.0000001;
            double twopi = Math.PI * 2;
            Point3D p1 = new Point3D();
            Point3D p2 = new Point3D();

            for (int i = 0; i < n; i++)
            {
                p1.X = p.points[i].X - q.X;
                p1.Y = p.points[i].Y - q.Y;
                p1.Z = p.points[i].Z - q.Z;
                p2.X = p.points[(i + 1) % n].X - q.X;
                p2.Y = p.points[(i + 1) % n].Y - q.Y;
                p2.Z = p.points[(i + 1) % n].Z - q.Z;

                m1 = Math.Sqrt(p1.X * p1.X + p1.Y * p1.Y + p1.Z * p1.Z);
                m2 = Math.Sqrt(p2.X * p2.X + p2.Y * p2.Y + p2.Z * p2.Z);
                if (m1 * m2 <= epsilon)
                    return (twopi);
                else
                    costheta = (p1.X * p2.X + p1.Y * p2.Y + p1.Z * p2.Z) / (m1 * m2);

                anglesum += Math.Acos(costheta);
            }
            return (anglesum);
        }

        public static Point Calc_X(Edge2D e1, Edge2D e2)
        {
            Point X = new Point(-1, -1);

            double a_1 = (e1.p2.Y - e1.p1.Y);
            double a_2 = (e2.p2.Y - e2.p1.Y);

            double b_1 = (e1.p1.X - e1.p2.X);
            double b_2 = (e2.p1.X - e2.p2.X);

            double c_1 = (e1.p1.X * (e1.p1.Y - e1.p2.Y) + e1.p1.Y * (e1.p2.X - e1.p1.X));
            double c_2 = (e2.p1.X * (e2.p1.Y - e2.p2.Y) + e2.p1.Y * (e2.p2.X - e2.p1.X));

            double D = a_1 * b_2 - b_1 * a_2;
            double D_X = (-c_1) * b_2 - b_1 * (-c_2);
            double D_Y = a_1 * (-c_2) - (-c_1) * a_2;

            X.X = (int)(D_X / D);
            X.Y = (int)(D_Y / D);

            return X;
        }

        public static bool Is_X_Edge(Edge2D e1, Edge2D e2)
        {
            float P1P2_X = e1.p2.X - e1.p1.X;
            float P3P4_X = e2.p2.X - e2.p1.X;

            float P1P3_X = e2.p1.X - e1.p1.X;
            float P1P4_X = e2.p2.X - e1.p1.X;

            float P3P1_X = e1.p1.X - e2.p1.X;
            float P3P2_X = e1.p2.X - e2.p1.X;

            float P1P2_Y = e1.p2.Y - e1.p1.Y;
            float P3P4_Y = e2.p2.Y - e2.p1.Y;

            float P1P3_Y = e2.p1.Y - e1.p1.Y;
            float P1P4_Y = e2.p2.Y - e1.p1.Y;

            float P3P1_Y = e1.p1.Y - e2.p1.Y;
            float P3P2_Y = e1.p2.Y - e2.p1.Y;

            float v1 = P3P4_X * P3P1_Y - P3P4_Y * P3P1_X;
            float v2 = P3P4_X * P3P2_Y - P3P4_Y * P3P2_X;
            float v3 = P1P2_X * P1P3_Y - P1P2_Y * P1P3_X;
            float v4 = P1P2_X * P1P4_Y - P1P2_Y * P1P4_X;

            int v1v2 = Math.Sign(v1) * Math.Sign(v2);
            int v3v4 = Math.Sign(v3) * Math.Sign(v4);

            if (v1v2 < 0 && v3v4 < 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        //Поиск пересечения
        public static bool find_intersection(Point3D point, Point3D dir, ref Point3D hit, ref Point3D N, ref Material material)
        {
            double dist = float.MaxValue;

            foreach (Mesh mesh in scene.mesh)
            {
                foreach (Polygon polygon in mesh.faces)
                {
                    Point3D p = new Point3D();
                    double dist_i = 0;
                    if (polygon.ray_intersection(point, dir, ref p, ref dist_i) && dist_i < dist)
                    {
                        dist = dist_i;
                        hit = point + dir * dist_i;
                        N = polygon.normal;
                        material = mesh.mat;
                    }
                }
            }

            foreach (Sphere sphere in scene.spheres)
            {
                double dist_i = 0;
                if (sphere.ray_intersection(point, dir, ref dist_i) && dist_i < dist)
                {
                    dist = dist_i;
                    hit = point + dir * dist_i;
                    N = normalize(hit - sphere.location);
                    material = sphere.mat;
                }
            }

            return dist < 1000000;
        }

        public static double Length(Point3D p1)
        {
            return Math.Sqrt(p1.X * p1.X + p1.Y * p1.Y + p1.Z * p1.Z);
        }

        public static double[,] LineRotate(Point3D p1, Point3D p2, double angle)
        {
            double l = (p2.X - p1.X) / Distance(p1, p2);
            double m = (p2.Y - p1.Y) / Distance(p1, p2);
            double n = (p2.Z - p1.Z) / Distance(p1, p2);

            return new double[4, 4]
                {{ l*l + Math.Cos(angle)*(1 - l*l), l*(1-Math.Cos(angle))*m + n*Math.Sin(angle), l*(1 - Math.Cos(angle))*n - m*Math.Sin(angle), 0 },
                { l*(1 - Math.Cos(angle))*m - n*Math.Sin(angle), m*m + Math.Cos(angle)*(1 - m*m), m*(1 - Math.Cos(angle))*n + l*Math.Sin(angle), 0 },
                { l*(1 - Math.Cos(angle))*n + m*Math.Sin(angle), m*(1 - Math.Cos(angle))*n - l*Math.Sin(angle), n*n + Math.Cos(angle)*(1 - n*n), 0 },
                { 0, 0, 0, 1}
                };
        }

        public static string LoadMeshWithPath(ref Mesh mesh, string path)
        {
            string ans = "object";

            string filename = path;
            string[] text = File.ReadAllLines(filename);
            mesh = new Mesh();
            int cnt = 0;
            int cntN = 0;
            List<Point3D> normals = new List<Point3D>();
            foreach (string x in text)
            {
                if (x.StartsWith("o "))
                {
                    ans = x.Remove(0, 2);
                }
                if (x.StartsWith("v "))
                {
                    string[] s = x.Remove(0, 2).Split(' ');
                    Point3D point = new Point3D();
                    point.X = Double.Parse(s[0], new CultureInfo("en-us"));
                    point.Y = Double.Parse(s[1], new CultureInfo("en-us"));
                    point.Z = Double.Parse(s[2], new CultureInfo("en-us"));
                    point.index = cnt;
                    mesh.points.Add(point);
                    cnt++;
                }

                if (x.StartsWith("f "))
                {
                    string[] s = x.Remove(0, 2).Split(' ');
                    List<Point3D> l = new List<Point3D>();
                    foreach (string s1 in s)
                    {
                        l.Add(mesh.points[int.Parse(s1.Split('/')[0]) - 1]);
                    }
                    Point3D norm = new Point3D(int.Parse(s[0].Split('/')[2]) - 1, 0, 0);
                    Polygon poly = new Polygon(l, norm);
                    mesh.faces.Add(poly);

                }
                if (x.StartsWith("vn "))
                {
                    string[] s = x.Remove(0, 3).Split(' ');
                    Point3D norm = new Point3D(
                        Convert.ToDouble(s[0], CultureInfo.InvariantCulture),
                        Convert.ToDouble(s[1], CultureInfo.InvariantCulture),
                        Convert.ToDouble(s[2], CultureInfo.InvariantCulture));
                    normals.Add(norm);
                    cntN++;
                }
            }
            foreach (Polygon pol in mesh.faces)
            {
                int ind = (int)pol.normal.X;
                pol.normal = normals[ind];
            }
            return ans;
        }

        public class Point3D : IComparable<Point3D>
        {
            public double X;
            public double Y;
            public double Z;
            public int index;
            public Point3D(double x, double y, double z, int ind)
            {
                X = x;
                Y = y;
                Z = z;
                index = ind;
            }
            public Point3D(double x, double y, double z)
            {
                X = x;
                Y = y;
                Z = z;
                index = 0;
            }
            public Point3D()
            {
                X = 0;
                Y = 0;
                Z = 0;
                index = 0;
            }
            public Point3D(Point3D p)
            {
                X = p.X;
                Y = p.Y;
                Z = p.Z;
                index = p.index;
            }
            public int CompareTo(Point3D that)
            {
                if (X < that.X) return -1;
                if (Y < that.Y) return -1;
                if (X == that.X && Y == that.Y) return 0;
                return 1;
            }
            public static double operator *(Point3D p1, Point3D p2)
            {
                return dotProduct(p1, p2);
            }
            public static Point3D operator +(Point3D p1, Point3D p2)
            {
                return new Point3D(p1.X + p2.X, p1.Y + p2.Y, p1.Z + p2.Z);
            }
            public static Point3D operator -(Point3D p1, Point3D p2)
            {
                return new Point3D(p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z);
            }
            public static Point3D operator *(Point3D p1, double d)
            {
                return new Point3D(p1.X * d, p1.Y * d, p1.Z * d);
            }
        }

        public class Edge
        {
            public Point3D p1;
            public Point3D p2;
            public Edge(Point3D pp1, Point3D pp2)
            {
                p1 = pp1;
                p2 = pp2;
            }
            public Edge()
            {
                p1 = new Point3D();
                p2 = new Point3D();
            }
        }

        public class Edge2D
        {
            public PointF p1;
            public PointF p2;
            public Edge2D(PointF pp1, PointF pp2)
            {
                p1 = pp1;
                p2 = pp2;
            }
            public Edge2D()
            {
                p1 = new PointF();
                p2 = new PointF();
            }
        }

        public class Polygon
        {
            public List<Point3D> points;
            public Point3D normal;
            public Polygon()
            {
                points = new List<Point3D>();
            }
            public Polygon(List<Point3D> l, bool clockwise = true)
            {
                points = new List<Point3D>();
                foreach (Point3D p in l)
                {
                    Point3D t = new Point3D(p);
                    points.Add(t);
                }
                normal = CreateNormal(this, clockwise);
            }
            public Polygon(List<Point3D> l, Point3D norm)
            {
                points = new List<Point3D>();
                foreach (Point3D p in l)
                {
                    Point3D t = new Point3D(p);
                    points.Add(t);
                }
                normal = new Point3D(norm);
            }

            //Было ли пересечение с плоскостью
            public bool ray_intersection(Point3D rayL, Point3D rayV, ref Point3D P, ref double dist)
            {
                Point3D Norm = normal;

                double d = PointPlaneDistance(this, new Point3D(0, 0, 0), false);
                double dot1 = rayL * Norm;
                double dot2 = rayV * Norm;

                if (Math.Abs(dot2) < 0.000001)
                {
                    return false;
                }

                double t = (-d - dot1) / dot2;
                dist = t;

                if (t < 0) return false;

                P = rayL + rayV * t;

                double twopi = Math.PI * 2;
                double anglesum = CalcAngleSum(P, this);

                if (Math.Abs(anglesum - twopi) > 0.0001)
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }
        }

        public class Polygon2D
        {
            public List<PointF> points;
            public Polygon2D()
            {
                points = new List<PointF>();
            }
            public Polygon2D(List<PointF> l)
            {
                points = new List<PointF>();
                foreach (PointF p in l)
                {
                    PointF t = new PointF(p.X, p.Y);
                    points.Add(t);
                }
            }

            public bool IsPointInside(PointF p)
            {
                Edge2D ray = new Edge2D(p, new PointF(p.X + 1000, p.Y));
                if (EdgeIntersectWithPoly(ray) % 2 == 0)
                    return false;
                else return true;

            }
            private int EdgeIntersectWithPoly(Edge2D e)
            {
                var arr_Point = points.ToArray();
                int intersect_counter = 0;
                for (int i = 0; i < arr_Point.Length; i++)
                {
                    Edge2D edge;
                    PointF p = new PointF();
                    if (i == arr_Point.Length - 1)
                        edge = new Edge2D(arr_Point[i], arr_Point[0]);
                    else
                        edge = new Edge2D(arr_Point[i], arr_Point[i + 1]);

                    if (CheckEdgesForIntersection(edge, e, ref p))
                        intersect_counter++;
                }
                return intersect_counter;
            }

            static public double MultVectors(PointF v1, PointF v2)
            {
                return v1.X * v2.Y - v1.Y * v2.X;
            }

            static public bool CheckEdgesForIntersection(Edge2D e1, Edge2D e2, ref PointF res)
            {
                float x1 = e1.p1.X;
                float y1 = e1.p1.Y;

                float x2 = e1.p2.X;
                float y2 = e1.p2.Y;

                float x3 = e2.p1.X;
                float y3 = e2.p1.Y;

                float x4 = e2.p2.X;
                float y4 = e2.p2.Y;

                PointF v_e2se2e = new PointF(x4 - x3, y4 - y3);
                PointF v_e2se1s = new PointF(x1 - x3, y1 - y3);
                PointF v_e2se1e = new PointF(x2 - x3, y2 - y3);
                PointF v_e1se1e = new PointF(x2 - x1, y2 - y1);
                PointF v_e1se2s = new PointF(x3 - x1, y3 - y1);
                PointF v_e1se2e = new PointF(x4 - x1, y4 - y1);

                double v1 = MultVectors(v_e2se2e, v_e2se1s);
                double v2 = MultVectors(v_e2se2e, v_e2se1e);
                double v3 = MultVectors(v_e1se1e, v_e1se2s);
                double v4 = MultVectors(v_e1se1e, v_e1se2e);

                double mult1 = v1 * v2;
                double mult2 = v3 * v4;

                if (mult1 < 0 && mult2 < 0)
                {
                    double a1 = y2 - y1;
                    double b1 = x1 - x2;
                    double c1 = x1 * (y1 - y2) + y1 * (x2 - x1);

                    double a2 = y4 - y3;
                    double b2 = x3 - x4;
                    double c2 = x3 * (y3 - y4) + y3 * (x4 - x3);
                    double det = a1 * b2 - a2 * b1;
                    double detx = c2 * b1 - c1 * b2;
                    double dety = c1 * a2 - a1 * c2;
                    res = new PointF((int)(detx / det), (int)(dety / det));
                    return true;
                }

                return false;

            }
        }

        public class Mesh
        {
            public List<Point3D> points;
            public SortedDictionary<int, List<int>> connections;
            public List<Edge> edges;
            public List<Polygon> faces;
            public Material mat;

            public Mesh()
            {
                points = new List<Point3D>();
                connections = new SortedDictionary<int, List<int>>();
                edges = new List<Edge>();
                faces = new List<Polygon>();
                mat = new Material();
            }

            public Mesh(List<Point3D> l, SortedDictionary<int, List<int>> sd, List<Edge> le, List<Polygon> lf)
            {
                points = new List<Point3D>();
                connections = new SortedDictionary<int, List<int>>();
                edges = new List<Edge>();
                faces = new List<Polygon>();
                mat = new Material();

                foreach (Point3D p in l)
                {
                    Point3D p3D = new Point3D(p);
                    points.Add(p3D);
                    List<int> temp = new List<int>();
                    if (sd.ContainsKey(p.index))
                        foreach (int pp in sd[p.index])
                        {
                            temp.Add(pp);
                        }
                    connections.Add(p.index, temp);
                }

                if (le.Count() == 0)
                {
                    int countP = points.Count();
                    bool[,] flags = new bool[countP, countP];

                    foreach (Point3D p1 in points)
                    {
                        int p1ind = p1.index;
                        foreach (int p2ind in connections[p1ind])
                        {
                            if (!flags[p1ind, p2ind])
                            {
                                flags[p1ind, p2ind] = true;
                                flags[p2ind, p1ind] = true;
                                Point3D t1 = new Point3D(p1);
                                Point3D t2 = new Point3D(points[p2ind]);
                                edges.Add(new Edge(t1, t2));
                            }
                        }
                    }
                }
                else
                {
                    foreach (Edge e in le)
                    {
                        Point3D t1 = new Point3D(e.p1);
                        Point3D t2 = new Point3D(e.p2);
                        edges.Add(new Edge(t1, t2));
                    }
                }
                if (lf.Count != 0)
                {
                    foreach (Polygon p in lf)
                    {
                        faces.Add(new Polygon(p.points));
                    }
                }
            }

            public Mesh(Mesh oldM)
            {
                var l = oldM.points;
                var sd = oldM.connections;
                var le = oldM.edges;
                var lf = oldM.faces;
                points = new List<Point3D>();
                connections = new SortedDictionary<int, List<int>>();
                edges = new List<Edge>();
                faces = new List<Polygon>();
                mat = new Material();
                foreach (Point3D p in l)
                {
                    Point3D p3D = new Point3D(p);
                    points.Add(p3D);
                    List<int> temp = new List<int>();
                    if (sd.ContainsKey(p.index))
                        foreach (int pp in sd[p.index])
                        {
                            temp.Add(pp);
                        }
                    connections.Add(p.index, temp);
                }
                if (le.Count() == 0)
                {
                    int countP = points.Count();
                    bool[,] flags = new bool[countP, countP];
                    foreach (Point3D p1 in points)
                    {
                        int p1ind = p1.index;
                        foreach (int p2ind in connections[p1ind])
                        {
                            if (!flags[p1ind, p2ind])
                            {
                                flags[p1ind, p2ind] = true;
                                flags[p2ind, p1ind] = true;
                                Point3D t1 = new Point3D(p1);
                                Point3D t2 = new Point3D(points[p2ind]);
                                edges.Add(new Edge(t1, t2));
                            }
                        }
                    }
                }
                else
                {
                    foreach (Edge e in le)
                    {
                        Point3D t1 = new Point3D(e.p1);
                        Point3D t2 = new Point3D(e.p2);
                        edges.Add(new Edge(t1, t2));
                    }
                }
                if (lf.Count != 0)
                {
                    foreach (Polygon p in lf)
                    {
                        faces.Add(new Polygon(p.points));
                    }
                }
            }

            public void Clear()
            {
                points.Clear();
                connections.Clear();
                edges.Clear();
                faces.Clear();
            }
        }

        //Отдельный класс для параметрисечкой сферы
        public class Sphere
        {
            public Point3D location;
            public double radius;
            public Material mat;
            public Sphere(Point3D l, double r)
            {
                location = new Point3D(l);
                radius = r;
                mat = new Material();
            }
            public Sphere()
            {
                location = new Point3D();
                radius = 0;
                mat = new Material();
            }
            //Было ли пересечение со сферой
            public bool ray_intersection(Point3D rayL, Point3D rayV, ref double dist)
            {
                Point3D L = new Point3D(location.X - rayL.X, location.Y - rayL.Y, location.Z - rayL.Z);
                double t = dotProduct(L, rayV);
                double d2 = dotProduct(L, L) - t * t;

                if (d2 > radius * radius)
                {
                    return false;
                }

                double thc = Math.Sqrt(radius * radius - d2);
                dist = t - thc;
                double t1 = t + thc;

                if (dist < 0)
                {
                    dist = t1;
                }

                if (dist < 0)
                {
                    return false;
                }
                else
                {
                    return true;
                }
            }
        }

        //Класс для света
        public class Light
        {
            public Point3D location;
            public double intensity;
            public Color color;

            public Light()
            {
                location = new Point3D();
                intensity = 0;
                color = Color.White;
            }
            public Light(Point3D p, double i, Color c)
            {
                location = new Point3D(p);
                intensity = i;
                color = c;
            }
        }

        //Класс материалов
        public class Material
        {
            public Color diffus;
            public double specular;
            public double refraction_index;
            public List<double> param;

            public Material(Color c, double spec, List<double> alb, double r)
            {
                diffus = c;
                specular = spec;
                param = new List<double> { alb[0], alb[1], alb[2], alb[3] };
                refraction_index = r;
            }
            public Material()
            {
                diffus = Color.Gray;
                specular = 0;
                param = new List<double> { 1, 0, 0 };
                refraction_index = 1;
            }
        }

        //Класс для сцены
        public class Scene
        {
            public List<Sphere> spheres;
            public List<Mesh> mesh;
            public List<Light> lights;

            public Scene()
            {
                spheres = new List<Sphere>();
                mesh = new List<Mesh>();
                lights = new List<Light>();
            }
        }

    }
}
