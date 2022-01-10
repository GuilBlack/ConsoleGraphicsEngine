//#include "Engine3D.h"
#include "vendor/olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <algorithm>

struct vec3d
{
    float x, y, z;
};

struct triangle
{
    vec3d p[3];

	wchar_t sym;
	short col;
};

struct mesh
{
    std::vector<triangle> triangles;

	bool LoadFromObjectFile(std::string sFileName)
	{
		std::ifstream f(sFileName);
		if (!f.is_open())
			return false;

		std::vector<vec3d> vertices;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				vec3d v;
				s >> junk >> v.x >> v.y >> v.z;
				vertices.push_back(v);
			}

			if (line[0] == 'f')
			{
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				triangles.push_back({ vertices[f[0] - 1], vertices[f[1] - 1], vertices[f[2] - 1] });
			}
		}

		return true;
	}
};

struct mat4x4f
{
    float mat[4][4] = { 0 };
};

class Engine3D : public olcConsoleGameEngine
{
public:
    Engine3D()
    {
        m_sAppName = L"3D Graphics Engine";
    }

private:
    mesh m_MeshCube;
    mat4x4f m_MatProj;

	vec3d m_vCamera = { 0, 0, 0 }; //simplified camera

	float m_fTheta;

    void MultiplyMatrixVector(const vec3d& in, vec3d &out, const mat4x4f& mat)
    {
        out.x = in.x * mat.mat[0][0] + in.y * mat.mat[1][0] + in.z * mat.mat[2][0] + mat.mat[3][0];
        out.y = in.x * mat.mat[0][1] + in.y * mat.mat[1][1] + in.z * mat.mat[2][1] + mat.mat[3][1];
        out.z = in.x * mat.mat[0][2] + in.y * mat.mat[1][2] + in.z * mat.mat[2][2] + mat.mat[3][2];
        float w = in.x * mat.mat[0][3] + in.y * mat.mat[1][3] + in.z * mat.mat[2][3] + mat.mat[3][3];

        if (w != 0.0f)
        {
            out.x /= w; out.y /= w; out.z /= w;
        }
    }

	// Taken From Command Line Webcam Video
	CHAR_INFO GetColor(float lum)
	{
		short bg_col, fg_col;
		wchar_t sym;
		int pixel_bw = (int)(13.0f * lum);
		switch (pixel_bw)
		{
		case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

		case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
		case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
		case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

		case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
		case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
		case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

		case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
		case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
		case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
		case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
		default:
			bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
		}

		CHAR_INFO c;
		c.Attributes = bg_col | fg_col;
		c.Char.UnicodeChar = sym;
		return c;
	}

public:
    bool OnUserCreate() override
	{
		//m_MeshCube.triangles = {

		//	// SOUTH
		//	{ 0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f },
		//	{ 0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		//	// EAST                                                      
		//	{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f },
		//	{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f },

		//	// NORTH                                                     
		//	{ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f },
		//	{ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f },

		//	// WEST                                                      
		//	{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f },
		//	{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f },

		//	// TOP                                                       
		//	{ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f },
		//	{ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f },

		//	// BOTTOM                                                    
		//	{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f },
		//	{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		//};

		m_MeshCube.LoadFromObjectFile("res/obj/LowPolyTree.obj");

		// Projection Matrix
		float fNear = 0.1f;
		float fFar = 1000.0f;
		float fFov = 90.0f;
		float fAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
		float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159f); //in radian and not in degrees

		m_MatProj.mat[0][0] = fAspectRatio * fFovRad;
		m_MatProj.mat[1][1] = fFovRad;
		m_MatProj.mat[2][2] = fFar / (fFar - fNear);
		m_MatProj.mat[3][2] = (-fFar * fNear) / (fFar - fNear);
		m_MatProj.mat[2][3] = 1.0f;
		m_MatProj.mat[3][3] = 0.0f;

		return true;
	}

    bool OnUserUpdate(float fElapseTime) override
	{
		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

		mat4x4f matRotZ, matRotX;
		m_fTheta += 1.0f * fElapseTime;

		//Rotation Z
		matRotZ.mat[0][0] = cosf(m_fTheta);
		matRotZ.mat[0][1] = sinf(m_fTheta);
		matRotZ.mat[1][0] = -sinf(m_fTheta);
		matRotZ.mat[1][1] = cosf(m_fTheta);
		matRotZ.mat[2][2] = 1;
		matRotZ.mat[3][3] = 1;

		//Rotation X
		matRotX.mat[0][0] = 1;
		matRotX.mat[1][1] = cosf(m_fTheta * 0.5f);
		matRotX.mat[1][2] = sinf(m_fTheta * 0.5f);
		matRotX.mat[2][1] = -sinf(m_fTheta * 0.5f);
		matRotX.mat[2][2] = cosf(m_fTheta * 0.5f);
		matRotX.mat[3][3] = 1;

		std::vector<triangle> trianglesToRaster;

		//Draw Triangles
		for (const auto &tri : m_MeshCube.triangles)
		{
			//rotates the cube
			triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;

			MultiplyMatrixVector(tri.p[0], triRotatedZ.p[0], matRotZ);
			MultiplyMatrixVector(tri.p[1], triRotatedZ.p[1], matRotZ);
			MultiplyMatrixVector(tri.p[2], triRotatedZ.p[2], matRotZ);

			MultiplyMatrixVector(triRotatedZ.p[0], triRotatedZX.p[0], matRotX);
			MultiplyMatrixVector(triRotatedZ.p[1], triRotatedZX.p[1], matRotX);
			MultiplyMatrixVector(triRotatedZ.p[2], triRotatedZX.p[2], matRotX);

			triTranslated = triRotatedZX;

			//offset it into our screen so that we can see
			triTranslated.p[0].z = triRotatedZX.p[0].z + 10.0f;
			triTranslated.p[1].z = triRotatedZX.p[1].z + 10.0f;
			triTranslated.p[2].z = triRotatedZX.p[2].z + 10.0f;

			vec3d normal, vecA, vecB; //normal is the normal vector to the plane made with the 2 vectors a and b
			vecA.x = triTranslated.p[1].x - triTranslated.p[0].x;
			vecA.y = triTranslated.p[1].y - triTranslated.p[0].y;
			vecA.z = triTranslated.p[1].z - triTranslated.p[0].z;

			vecB.x = triTranslated.p[2].x - triTranslated.p[0].x;
			vecB.y = triTranslated.p[2].y - triTranslated.p[0].y;
			vecB.z = triTranslated.p[2].z - triTranslated.p[0].z;

			normal.x = vecA.y * vecB.z - vecA.z * vecB.y;
			normal.y = vecA.z * vecB.x - vecA.x * vecB.z;
			normal.z = vecA.x * vecB.y - vecA.y * vecB.x;

			float normalLength = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z); // returns the length of the normal
			normal.x /= normalLength; normal.y /= normalLength; normal.z /= normalLength; //normalize the normal

			//vector between the surface normal and the camera
			vec3d projLine = {
				triTranslated.p[0].x - m_vCamera.x,
				triTranslated.p[0].y - m_vCamera.y,
				triTranslated.p[0].z - m_vCamera.z,
			};

			//negative only if the 2 vectors are facing different directions
			float dotProdNormProjLine = normal.x*projLine.x + normal.y*projLine.y + normal.z*projLine.z;

			//draw triangle only if the normal is facing the cam so if its z is negative
			if (dotProdNormProjLine < 0)
			{
				//Lighting
				vec3d directionalLight = { 0.0f, 0.0f, -1.0f }; //light facing the the cube
				float lightLength = sqrtf(directionalLight.x*directionalLight.x + directionalLight.y*directionalLight.y + 
					directionalLight.z*directionalLight.z);
				directionalLight.x /= lightLength; directionalLight.y /= lightLength; directionalLight.z /= lightLength;

				float dotProdNormLight = normal.x*directionalLight.x + normal.y*directionalLight.y + normal.z*directionalLight.z;

				CHAR_INFO c = GetColor(dotProdNormLight);
				triTranslated.col = c.Attributes;
				triTranslated.sym = c.Char.UnicodeChar;

				//projects the triangles from 3d to 2d
				MultiplyMatrixVector(triTranslated.p[0], triProjected.p[0], m_MatProj);
				MultiplyMatrixVector(triTranslated.p[1], triProjected.p[1], m_MatProj);
				MultiplyMatrixVector(triTranslated.p[2], triProjected.p[2], m_MatProj);
				triProjected.col = triTranslated.col;
				triProjected.sym = triTranslated.sym;

				//Scale into view

				//will make range -1 / 1 to 0 / 2
				triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
				triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
				triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;

				//rescale it between 0 and 1 then adjusting it to screen sizes
				triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
				triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
				triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
				triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

				trianglesToRaster.push_back(triProjected);
			}

		}

		//based on the painter algorithm to draw the furthest triangles first
		std::sort(trianglesToRaster.begin(), trianglesToRaster.end(), 
			[](triangle &t1, triangle &t2)
			{
				//taking the z of the mid point of the faces and comparing them
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
				return z1 > z2;
			});

		for (auto& triProjected : trianglesToRaster)
		{
			FillTriangle((int)triProjected.p[0].x, (int)triProjected.p[0].y,
				(int)triProjected.p[1].x, (int)triProjected.p[1].y,
				(int)triProjected.p[2].x, (int)triProjected.p[2].y,
				triProjected.sym, triProjected.col);

			DrawTriangle((int)triProjected.p[0].x, (int)triProjected.p[0].y,
				(int)triProjected.p[1].x, (int)triProjected.p[1].y,
				(int)triProjected.p[2].x, (int)triProjected.p[2].y,
				PIXEL_SOLID, FG_WHITE);
		}

		return true;
	}
};

int main()
{
    Engine3D engine;
	if (engine.ConstructConsole(256, 156, 4, 4))
		engine.Start();
	else
		std::cout << "[ERROR]: Console didn't start" << std::endl;

	return 0;
}
