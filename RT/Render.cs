using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Drawing;
using System.Windows.Forms;

using static RT.Geometry;
using static RT.Render;

namespace RT
{
    class Render
    {
        public static Point3D Camera;

        public static Color Ray_Cast(Point3D point, Point3D dir, int depth = 0)
        {
            //Точка пересечения луча и меша
            Point3D hit = new Point3D();
            //нормаль в точке
            Point3D N = new Point3D();

            //Материал в точке пересечения
            Material mat = new Material();
            //Дистанция до обьекта (до ближайщего обьекта)
            double dist = double.MaxValue;

            //Было ли пересечение луча света с обьектами сцены
            bool intersect = false;

            //Ищем пересечение со всеми мешами
            foreach (Mesh mesh in scene.mesh)
            {
                foreach (Polygon pol in mesh.faces)
                {

                    double a = 0;
                    Point3D P = new Point3D();
                    if (depth < 4 && pol.ray_intersection(point, dir, ref P, ref a))
                    {
                        if (a < dist)
                        {
                            intersect = true;

                            dist = a;
                            hit = P;
                            N = pol.normal;
                            mat = mesh.mat;
                        }
                    }
                }
            }

            //Ищем пересечение со всеми сферами
            foreach (Sphere sphere in scene.spheres)
            {
                double a = 0;
                if (depth < 4 && sphere.ray_intersection(point, dir, ref a))
                {
                    if (a < dist)
                    {
                        intersect = true;

                        hit = point + dir * a;
                        N = normalize(hit - sphere.location);
                        dist = a;

                        mat = sphere.mat;
                    }
                }
            }

            if (intersect)
            {
                //Интенсивность диффузного и спекуляр освещения
                double diffuse_light_intensity = 0;
                double specular_light_intensity = 0.0;

                //Цвета отражения/приломления
                Color reflect_color = new Color();
                Color refract_color = new Color();

                //Каждый источник света обрабатывается отдельно, а потом результаты смешиваются
                foreach (Light light in scene.lights)
                {
                    //Для теней
                    Point3D lightVecFull = light.location - hit;
                    Point3D lightVecN = normalize(lightVecFull);
                    double light_distance = Length(lightVecFull);

                    Point3D shadow_point = lightVecN * N < 0 ? hit - N * 0.001 : hit + N * 0.001;
                    Point3D shadow_hit = new Point3D();
                    Point3D shadow_N = new Point3D();
                    Material mat_s = new Material();

                    //Если материал стеклянный, тень на обьект не отбрасывается
                    if (mat.param[3] == 0)
                    {
                        if (find_intersection(shadow_point, lightVecN, ref shadow_hit, ref shadow_N, ref mat_s) && Length(shadow_hit - shadow_point) < light_distance)
                        {
                            continue;
                        }
                    }

                    //Освещение
                    diffuse_light_intensity += light.intensity * Math.Max(0.0, (lightVecN * N));
                    specular_light_intensity += Math.Pow(Math.Max(0.0, reflect(lightVecN, N) * dir), mat.specular) * light.intensity;

                    //Просчитываем, только если материал зеркальный
                    if (mat.param[2] != 0)
                    {
                        Point3D reflect_dir = normalize(reflect(dir, N));
                        Point3D reflect_orig = reflect_dir * N < 0 ? hit - N * 0.001 : hit + N * 0.001;

                        reflect_color = Ray_Cast(reflect_orig, reflect_dir, depth + 1);
                    }

                    //Просчитываем только если материал прозрачный
                    if (mat.param[3] != 0)
                    {
                        Point3D refract_dir = normalize(refract(dir, N, mat.refraction_index));
                        Point3D refract_orig = refract_dir * N < 0 ? hit - N * 0.001 : hit + N * 0.001;

                        refract_color = Ray_Cast(refract_orig, refract_dir, depth + 1);
                    }
                }

                double diffuseR = 0;
                double diffuseG = 0;
                double diffuseB = 0;

                double specularR = 0;
                double specularG = 0;
                double specularB = 0;

                double reflectR = 0;
                double reflectG = 0;
                double reflectB = 0;

                double refractR = 0;
                double refractG = 0;
                double refractB = 0;

                //Считаем все параметры, только если их альбедо не нулевой
                if (mat.param[0] != 0)
                {
                    diffuseR = mat.diffus.R * diffuse_light_intensity;
                    diffuseG = mat.diffus.G * diffuse_light_intensity;
                    diffuseB = mat.diffus.B * diffuse_light_intensity;
                }

                if (mat.param[1] != 0)
                {
                    specularR = mat.diffus.R * diffuse_light_intensity * mat.param[0] * specular_light_intensity * mat.param[1];
                    specularG = mat.diffus.G * diffuse_light_intensity * mat.param[0] * specular_light_intensity * mat.param[1];
                    specularB = mat.diffus.B * diffuse_light_intensity * mat.param[0] * specular_light_intensity * mat.param[1];
                }

                if (mat.param[2] != 0)
                {
                    reflectR = reflect_color.R * mat.param[2];
                    reflectG = reflect_color.G * mat.param[2];
                    reflectB = reflect_color.B * mat.param[2];
                }

                if (mat.param[3] != 0)
                {
                    refractR = refract_color.R * mat.param[3];
                    refractG = refract_color.G * mat.param[3];
                    refractB = refract_color.B * mat.param[3];
                }

                int r = (int)(diffuseR + specularR + reflectR + refractR);
                int g = (int)(diffuseG + specularG + reflectG + refractG);
                int b = (int)(diffuseB + specularB + reflectB + refractB);

                //Защита от ошибок вычислений
                if (r > 255)
                {
                    r = 255;
                }

                if (g > 255)
                {
                    g = 255;
                }

                if (b > 255)
                {
                    b = 255;
                }

                //Чтобы были более яркие блики
                if (specular_light_intensity > 1)
                {
                    r = 255;
                    g = 255;
                    b = 255;
                }

                return Color.FromArgb(r, g, b);
            }

            //Цвет фона
            return Color.FromArgb(0, 0, 0);
        }

        public static void Render_Func(Bitmap bm, PictureBox pb)
        {
            Graphics g = Graphics.FromImage(bm);
            g.Clear(Color.Transparent);

            //Выпускаем луч из каждого пикселя. Обратная трассировка лучей.
            for (int i = 0; i < bm.Width; i++)
            {
                for (int j = 0; j < bm.Height; j++)
                {
                    int x = i - bm.Width / 2;
                    int y = j - bm.Height / 2;
                    int z = (int)(-Camera.Z);

                    double distance = Distance(new Point3D(x, y, z), new Point3D(0, 0, 0));

                    bm.SetPixel(i, j, Ray_Cast(Camera, new Point3D(x / distance, y / distance, z / distance)));
                }
            }

            pb.Refresh();
        }
    }
}
