using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

using static RT.Geometry;
using static RT.Render;
using static RT.Athens;

namespace RT
{
    public partial class Form1 : Form
    {
        Bitmap pic;

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            Camera = new Point3D(0, 0, -500, 0);
            pic = new Bitmap(RenderTarget.Width, RenderTarget.Height);

            RenderTarget.Image = pic;
            scene = new Scene();

            LoadScene();
            LoadRoom();

            Render_Func(pic, RenderTarget);
        }

        public static void LoadScene()
        {
            //var octa_1 = new Mesh();
            //LoadMeshWithPath(ref octa_1, "../../models/octa.obj");
            //octa_1.mat = new Material();
            //octa_1.mat.diffus = Color.White;
            //octa_1.mat.param = new List<double> { 1.0, 0.0, 0.0, 0.0 };
            //AtheneTransform(ref octa_1, AtheneScale(0.75, 0.75, 0.75));
            //AtheneTransform(ref octa_1, AtheneRotate(-30.0 / 360.0 * (2.0 * Math.PI), 'y'), true);
            //AtheneTransform(ref octa_1, AtheneMove(0, -50, 50));
            //scene.mesh.Add(octa_1);


            //var octa_2 = new Mesh();
            //LoadMeshWithPath(ref octa_2, "../../models/octa.obj");
            //octa_2.mat = new Material();
            //octa_2.mat.diffus = Color.White;
            //octa_2.mat.refraction_index = 1.3;
            //octa_2.mat.param = new List<double> { 0.0, 0.0, 0.0, 0.9 };
            //AtheneTransform(ref octa_2, AtheneScale(0.4, 0.4, 0.4));
            //AtheneTransform(ref octa_2, AtheneRotate(30.0 / 360.0 * (2.0 * Math.PI), 'y'), true);
            //AtheneTransform(ref octa_2, AtheneMove(-30, 80, 30));
            //scene.mesh.Add(octa_2);

            var cube_1 = new Mesh();
            LoadMeshWithPath(ref cube_1, "../../models/cube.obj");
            cube_1.mat = new Material();
            cube_1.mat.diffus = Color.White;
            cube_1.mat.refraction_index = 1.3;
            cube_1.mat.param = new List<double> { 0.0, 0.0, 0.0, 0.5 };
            AtheneTransform(ref cube_1, AtheneScale(0.3, 0.3, 0.3));
            AtheneTransform(ref cube_1, AtheneRotate(-5.0 / 360.0 * (2.0 * Math.PI), 'y'), true);
            AtheneTransform(ref cube_1, AtheneRotate(0.0 / 360.0 * (2.0 * Math.PI), 'x'), true);
            AtheneTransform(ref cube_1, AtheneMove(-90, -130, -30));
            scene.mesh.Add(cube_1);

            var cube_2 = new Mesh();
            LoadMeshWithPath(ref cube_2, "../../models/cube.obj");
            cube_2.mat = new Material();
            cube_2.mat.diffus = Color.YellowGreen;
            cube_2.mat.param = new List<double> { 1.0, 0.0, 0.0, 0.0 };
            AtheneTransform(ref cube_2, AtheneScale(0.35, 0.35, 0.35));
            AtheneTransform(ref cube_2, AtheneRotate(30.0 / 360.0 * (2.0 * Math.PI), 'y'), true);
            AtheneTransform(ref cube_2, AtheneMove(100, 120, -120));
            scene.mesh.Add(cube_2);

            var sphere_1 = new Sphere();
            sphere_1.mat = new Material();
            sphere_1.mat.diffus = Color.Yellow;
            sphere_1.mat.specular = 50;
            sphere_1.mat.refraction_index = 1.3;
            sphere_1.mat.param = new List<double> { 1.0, 0.0, 0.0, 0.0 };
            sphere_1.location = new Point3D(-130, 30, 200);
            sphere_1.radius = 50;
            scene.spheres.Add(sphere_1);

            //var sphere_2 = new Sphere();
            //sphere_2.mat = new Material();
            //sphere_2.mat.diffus = Color.Red;
            //sphere_2.mat.specular = 1;
            //sphere_2.mat.param = new List<double> { 0.9, 0.1, 0.0, 0.0 };
            //sphere_2.location = new Point3D(-115, 120, -50);
            //sphere_2.radius = 30;
            //scene.spheres.Add(sphere_2);

            //var sphere_3 = new Sphere();
            //sphere_3.mat = new Material();
            //sphere_3.mat.diffus = Color.White;
            //sphere_3.mat.refraction_index = 1.3;
            //sphere_3.mat.param = new List<double> { 0.0, 0.0, 0.0, 0.9 };
            //sphere_3.location = new Point3D(90, 35, -40);
            //sphere_3.radius = 40;
            //scene.spheres.Add(sphere_3);

            //var sphere_4 = new Sphere();
            //sphere_4.mat = new Material();
            //sphere_4.mat.diffus = Color.White;
            //sphere_4.mat.param = new List<double> { 0.0, 0.0, 1.0, 0.0 };
            //sphere_4.location = new Point3D(85, -155, 100);
            //sphere_4.radius = 50;
            //scene.spheres.Add(sphere_4);

            var light_1 = new Light();
            light_1.location = new Point3D(150, 0, 0);
            light_1.intensity = 0.5;
            light_1.color = Color.White;
            scene.lights.Add(light_1);

            var light_2 = new Light();
            light_2.location = new Point3D(-150, 0, 0);
            light_2.intensity = 0.3;
            light_2.color = Color.White;
            scene.lights.Add(light_2);
        }

        private void LoadRoom()
        {
            var floor = new Mesh();
            LoadMeshWithPath(ref floor, "../../models/plane.obj");
            floor.mat = new Material();
            floor.mat.diffus = Color.Orange;
            floor.mat.param = new List<double> { 1.0, 0.0, 0.0, 0.0 };
            AtheneTransform(ref floor, AtheneScale(2.5, 2.5, 5.0));
            AtheneTransform(ref floor, AtheneRotate((180) / 360.0 * (2.0 * Math.PI), 'z'), true);
            AtheneTransform(ref floor, AtheneMove(0, -100 + 250, -1));
            scene.mesh.Add(floor);

            var ceiling = new Mesh();
            LoadMeshWithPath(ref ceiling, "../../models/plane.obj");
            ceiling.mat = new Material();
            ceiling.mat.diffus = Color.ForestGreen;
            ceiling.mat.param = new List<double> { 1.0, 0.0, 0.0, 0.0 };
            AtheneTransform(ref ceiling, AtheneScale(2.5, 2.5, 5.0));
            AtheneTransform(ref ceiling, AtheneMove(0, -100 - 250, -1));
            scene.mesh.Add(ceiling);

            var wall_1 = new Mesh();
            LoadMeshWithPath(ref wall_1, "../../models/plane.obj");
            wall_1.mat = new Material();
            wall_1.mat.diffus = Color.Yellow;
            wall_1.mat.param = new List<double> { 1.0, 0.0, 0.0, 0.0 };
            AtheneTransform(ref wall_1, AtheneScale(2.5, 2.5, 5.0));
            AtheneTransform(ref wall_1, AtheneRotate((90) / 360.0 * (2.0 * Math.PI), 'z'), true);
            AtheneTransform(ref wall_1, AtheneMove(-250, -100, -1));
            scene.mesh.Add(wall_1);

            var wall_2 = new Mesh();
            LoadMeshWithPath(ref wall_2, "../../models/plane.obj");
            wall_2.mat = new Material();
            wall_2.mat.diffus = Color.Red;
            wall_2.mat.refraction_index = 1.3;
            wall_2.mat.param = new List<double> { 0.0, 0.0, 1.0, 0.1 };
            AtheneTransform(ref wall_2, AtheneScale(2.5, 2.5, 2.5));
            AtheneTransform(ref wall_2, AtheneRotate((90) / 360.0 * (2.0 * Math.PI), 'x'), true);
            AtheneTransform(ref wall_2, AtheneMove(0, -100, -1 + 500));
            scene.mesh.Add(wall_2);

            var wall_3 = new Mesh();
            LoadMeshWithPath(ref wall_3, "../../models/plane.obj");
            wall_3.mat = new Material();
            wall_3.mat.diffus = Color.DarkOrchid;
            wall_3.mat.specular = 100;
            wall_3.mat.param = new List<double> { 1.0, 0.0, 0.0, 0.0 };
            AtheneTransform(ref wall_3, AtheneScale(2.5, 2.5, 5.0));
            AtheneTransform(ref wall_3, AtheneRotate((-90) / 360.0 * (2.0 * Math.PI), 'z'), true);
            AtheneTransform(ref wall_3, AtheneMove(250, -100, -1));
            scene.mesh.Add(wall_3);

            var wall_4 = new Mesh();
            LoadMeshWithPath(ref wall_4, "../../models/plane.obj");
            wall_4.mat = new Material();
            wall_4.mat.diffus = Color.White;
            wall_4.mat.param = new List<double> { 1.0, 0.0, 0.0, 0.0 };
            AtheneTransform(ref wall_4, AtheneScale(2.5, 2.5, 2.5));
            AtheneTransform(ref wall_4, AtheneRotate((-90) / 360.0 * (2.0 * Math.PI), 'x'), true);
            AtheneTransform(ref wall_4, AtheneMove(0, -100, -500 - 1));
            scene.mesh.Add(wall_4);
        }
    }
}
