/*
A 3D implementation of "A velocity correcting method for volume preserving viscoelastic fluids." [1], using the Shape Matching [2] velocity correction scheme. For now, it is using the SPH algorithm presented by Müller et al.[5], instead of IISPH [3]. Will most likely implement Divergence-free SPH [4] since it requires less memory than IISPH, as well as being more stable.

[1] Takahashi, Tetsuya, Issei Fujishiro, and Tomoyuki Nishita. "A velocity correcting method for volume preserving viscoelastic fluids." Proceedings of the Computer Graphics International. 2014.
[2] Müller, Matthias, et al. "Meshless deformations based on shape matching." ACM transactions on graphics (TOG) 24.3 (2005): 471-478.
[3] Ihmsen, Markus, et al. "Implicit incompressible SPH." IEEE Transactions on Visualization and Computer Graphics 20.3 (2014): 426-435.
[4] Bender, Jan, and Dan Koschier. "Divergence-free SPH for incompressible and viscous fluids." IEEE transactions on visualization and computer graphics 23.3 (2017): 1193-1206.
[5] Müller, Matthias, David Charypar, and Markus Gross. "Particle-based fluid simulation for interactive applications." Proceedings of the 2003 ACM SIGGRAPH/Eurographics symposium on Computer animation. Eurographics Association, 2003.

@author Octavio Navarro
@version 1.0
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <iostream>

#include <SPH_SM.h>

#if defined(_WIN32) || defined(_WIN64)
#  define WINDOWS_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#endif

// OpenGL Graphics includes
// #include <GL/glew.h>
#if defined (__APPLE__) || defined(MACOSX)
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/freeglut.h>
#endif

using namespace std;

// mouse controls
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -3.0;

bool keypressed = false, simulate = true;

float cubeFaces[24][3] = 
{
	{0.0,0.0,0.0}, {1.0,0.0,0.0}, {1.0,1.0,0.0}, {0.0,1.0,0.0},
	{0.0,0.0,0.0}, {1.0,0.0,0.0}, {1.0,0,1.0}, {0.0,0,1.0},
	{0.0,0.0,0.0}, {0,1.0,0.0}, {0,1.0,1.0}, {0.0,0,1.0},
	{0.0,1.0,0.0}, {1.0,1.0,0.0}, {1.0,1.0,1.0}, {0.0,1.0,1.0},
	{1.0,0.0,0.0}, {1.0,1.0,0.0}, {1.0,1.0,1.0}, {1.0,0.0,1.0},
	{0.0,0.0,1.0}, {0,1.0,1.0}, {1.0,1.0,1.0}, {1.0,0,1.0}
};

SPH_SM *sph;
int winX = 600;
int winY = 600;

void display_cube()
{
	glPushMatrix();
		glColor3f(1.0, 1.0, 1.0);
		glLineWidth(2.0);
		for(int i = 0; i< 6; i++)
		{
			glBegin(GL_LINE_LOOP);
			for(int j = 0; j < 4; j++)
				glVertex3f(cubeFaces[i*4+j][0], cubeFaces[i*4+j][1], cubeFaces[i*4+j][2]);
			glEnd();
		}
	glPopMatrix();
}

void display_points()
{
	glPushMatrix();
		Particle *p = sph->Get_Paticles();
		glColor3f(0.2f, 0.5f, 1.0f);
		glPointSize(2.0f);

		glBegin(GL_POINTS);
		for(int i=0; i<sph->Get_Particle_Number(); i++)
			glVertex3f(p[i].pos.x, p[i].pos.y, p[i].pos.z);
		glEnd();
	glPopMatrix();
}

void display (void)
{
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    glTranslatef(0.0, 0.0, translate_z);
    glRotatef(rotate_x, 1.0, 0.0, 0.0);
    glRotatef(rotate_y, 0.0, 1.0, 0.0);
	
	glPushMatrix();
		glLineWidth(10.0);
		glBegin(GL_LINES);
		glColor3f(0, 0, 1);
		glVertex3f(0, 0, 0);
		glVertex3f(1, 0, 0);

		glColor3f(1, 0, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 1, 0);

		glColor3f(0, 1, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 1);
		glEnd();
	glPopMatrix();

	display_cube();
	display_points();
	
	glutSwapBuffers();
}

void idle(void)
{
	if(simulate)
		sph->Animation();
	glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN) {
        mouse_buttons |= 1<<button;
    } else if (state == GLUT_UP) {
        mouse_buttons = 0;
    }

    mouse_old_x = x;
    mouse_old_y = y;
}

void motion(int x, int y)
{
    float dx, dy;
    dx = (float)(x - mouse_old_x);
    dy = (float)(y - mouse_old_y);

    if (mouse_buttons & 1) {
        rotate_x += dy * 0.2f;
        rotate_y += dx * 0.2f;
    } else if (mouse_buttons & 4) {
        translate_z += dy * 0.01f;
    }

    mouse_old_x = x;
    mouse_old_y = y;
}

void keys (unsigned char key, int x, int y)
{
	switch (key) {
		case 27:
			delete sph;
            exit(0);
			break;
		case 'q':
			cout << "Qm: " << sph->flip_quadratic() << endl;
			break;
		case 'v':
			cout << "Vc: " << sph->flip_volume() << endl;
			break;
		case 32:
			simulate = !simulate;
			cout << "Streaming: " << simulate << endl;
			break;
	}
}

void reshape (int w, int h)
{
	glViewport (0,0, w,h);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	
	gluPerspective (45, w*1.0/h*1.0, 0.01, 400);
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
}

void initGL ()
{
	const GLubyte* renderer;
	const GLubyte* version;
	const GLubyte* glslVersion;

	renderer = glGetString(GL_RENDERER); /* get renderer string */
	version = glGetString(GL_VERSION); /* version as a string */
	glslVersion = glGetString(GL_SHADING_LANGUAGE_VERSION);

	printf("Renderer: %s\n", renderer);
	printf("OpenGL version supported %s\n", version);
	printf("GLSL version supported %s\n", glslVersion);

	glEnable(GL_DEPTH_TEST);
}

void init(void)
{
	sph = new SPH_SM();

	std::vector<m3Vector> positions;
	m3Vector World_Size = m3Vector(1.0f, 1.0f, 1.0f);
	float kernel = 0.04f;

	for(float k = World_Size.z * 0.3f; k < World_Size.z * 0.7f; k += kernel * 0.6f)
	for(float j = World_Size.y * 0.0f; j < World_Size.y * 0.4f; j += kernel * 0.6f)
	for(float i = World_Size.x * 0.3f; i < World_Size.x * 0.7f; i += kernel * 0.6f)
		positions.push_back(m3Vector(i, j, k));

	sph->Init_Fluid(positions);
}

int main( int argc, const char **argv ) {

	srand((unsigned int)time(NULL));
	glutInit(&argc, (char**)argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize (winX, winY); 
	glutInitWindowPosition (100, 100);
	glutCreateWindow ("SPH SM 3D");
	glutReshapeFunc (reshape);
	glutKeyboardFunc (keys);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);

	initGL ();
	init();

	glutDisplayFunc(display); 
	glutIdleFunc (idle);
    glutMainLoop();

	return 0;
}
