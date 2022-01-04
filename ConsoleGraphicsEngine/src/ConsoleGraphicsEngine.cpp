//#include "Engine3D.h"
#include "vendor/olcConsoleGameEngine.h"

struct vec3d
{
    float x, y, z;
};

struct triangle
{
    vec3d p[3];
};

struct mesh
{
    std::vector<triangle> triangles;
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

public:
    bool OnUserCreate() override
	{
		m_MeshCube.triangles = {

			// SOUTH
			{ 0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f },
			{ 0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f },

			// EAST                                                      
			{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f },
			{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f },

			// NORTH                                                     
			{ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f },
			{ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f },

			// WEST                                                      
			{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f },
			{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f },

			// TOP                                                       
			{ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f },
			{ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f },

			// BOTTOM                                                    
			{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f },
			{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		};

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

		//Draw Triangles
		for (const auto &tri : m_MeshCube.triangles)
		{
			triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;

			MultiplyMatrixVector(tri.p[0], triRotatedZ.p[0], matRotZ);
			MultiplyMatrixVector(tri.p[1], triRotatedZ.p[1], matRotZ);
			MultiplyMatrixVector(tri.p[2], triRotatedZ.p[2], matRotZ);

			MultiplyMatrixVector(triRotatedZ.p[0], triRotatedZX.p[0], matRotX);
			MultiplyMatrixVector(triRotatedZ.p[1], triRotatedZX.p[1], matRotX);
			MultiplyMatrixVector(triRotatedZ.p[2], triRotatedZX.p[2], matRotX);

			triTranslated = triRotatedZX;
			triTranslated.p[0].z = triRotatedZX.p[0].z + 3.0f;
			triTranslated.p[1].z = triRotatedZX.p[1].z + 3.0f;
			triTranslated.p[2].z = triRotatedZX.p[2].z + 3.0f;

			MultiplyMatrixVector(triTranslated.p[0], triProjected.p[0], m_MatProj);
			MultiplyMatrixVector(triTranslated.p[1], triProjected.p[1], m_MatProj);
			MultiplyMatrixVector(triTranslated.p[2], triProjected.p[2], m_MatProj);

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
