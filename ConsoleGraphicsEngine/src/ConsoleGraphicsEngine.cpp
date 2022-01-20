//#include "Engine3D.h"
#include "vendor/olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <algorithm>

struct mat4x4f
{
    float mat[4][4] = { 0 };

	#pragma region mat4x4fUtilityMatrices

	static mat4x4f CreateRotationX(float fAngleRad)
	{
		mat4x4f matrix;
		matrix.mat[0][0] = 1.0f;
		matrix.mat[1][1] = cosf(fAngleRad);
		matrix.mat[1][2] = sinf(fAngleRad);
		matrix.mat[2][1] = -sinf(fAngleRad);
		matrix.mat[2][2] = cosf(fAngleRad);
		matrix.mat[3][3] = 1.0f;
		return matrix;
	}

	static mat4x4f CreateRotationY(float fAngleRad)
	{
		mat4x4f matrix;
		matrix.mat[0][0] = cosf(fAngleRad);
		matrix.mat[0][2] = sinf(fAngleRad);
		matrix.mat[2][0] = -sinf(fAngleRad);
		matrix.mat[1][1] = 1.0f;
		matrix.mat[2][2] = cosf(fAngleRad);
		matrix.mat[3][3] = 1.0f;
		return matrix;
	}

	static mat4x4f CreateRotationZ(float fAngleRad)
	{
		mat4x4f matrix;
		matrix.mat[0][0] = cosf(fAngleRad);
		matrix.mat[0][1] = sinf(fAngleRad);
		matrix.mat[1][0] = -sinf(fAngleRad);
		matrix.mat[1][1] = cosf(fAngleRad);
		matrix.mat[2][2] = 1.0f;
		matrix.mat[3][3] = 1.0f;
		return matrix;
	}

	static mat4x4f CreateIdentity()
	{
		mat4x4f matrix;
		matrix.mat[0][0] = 1.0f;
		matrix.mat[1][1] = 1.0f;
		matrix.mat[2][2] = 1.0f;
		matrix.mat[3][3] = 1.0f;
		return matrix;
	}

	static mat4x4f CreateTranslation(float x, float y, float z)
	{
		mat4x4f matrix;
		matrix.mat[0][0] = 1.0f;
		matrix.mat[1][1] = 1.0f;
		matrix.mat[2][2] = 1.0f;
		matrix.mat[3][3] = 1.0f;
		matrix.mat[3][0] = x;
		matrix.mat[3][1] = y;
		matrix.mat[3][2] = z;
		return matrix;
	}

	static mat4x4f CreateProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
	{
		float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
		mat4x4f matrix;
		matrix.mat[0][0] = fAspectRatio * fFovRad;
		matrix.mat[1][1] = fFovRad;
		matrix.mat[2][2] = fFar / (fFar - fNear);
		matrix.mat[3][2] = (-fFar * fNear) / (fFar - fNear);
		matrix.mat[2][3] = 1.0f;
		matrix.mat[3][3] = 0.0f;
		return matrix;
	}

	#pragma endregion


	#pragma region mat4x4fOperators

	mat4x4f& operator*=(const mat4x4f& rhs)
	{
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				this->mat[r][c] = this->mat[r][0] * rhs.mat[0][c] + this->mat[r][1] * rhs.mat[1][c] + this->mat[r][2] * rhs.mat[2][c] +
				this->mat[r][3] * rhs.mat[3][c];

		return *this;
	}

	mat4x4f operator*(const mat4x4f& rhs) const
	{
		mat4x4f matrix;
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				matrix.mat[r][c] = this->mat[r][0] * rhs.mat[0][c] + this->mat[r][1] * rhs.mat[1][c] + this->mat[r][2] * rhs.mat[2][c] +
				this->mat[r][3] * rhs.mat[3][c];

		return matrix;
	}

	#pragma endregion

};

struct vec3d
{
    float x, y, z, w;

	vec3d()
	{
		x = y = z = 0;
		w = 1;
	}

	vec3d(float a, float b, float c)
	{
		x = a; y = b; z = c; w = 1;
	}

	float GetLength()
	{
		return sqrtf(x*x + y*y + z*z);
	}

	vec3d Normalize()
	{
		float l = GetLength();
		return { x / l, y / l, z / l };
	}

	float DotProduct(vec3d& b)
	{
		return this->x * b.x + this->y * b.y + this->z * b.z;
	}

	vec3d CrossProduct(vec3d& b)
	{
		vec3d v;
		v.x = this->y * b.z - this->z * b.y;
		v.y = this->z * b.x - this->x * b.z;
		v.z = this->x * b.y - this->y * b.x;
		return v;
	}

	#pragma region vec3dOperators

	vec3d& operator+=(const vec3d& rhs)
	{
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	}

	vec3d& operator-=(const vec3d& rhs)
	{
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
		return *this;
	}

	vec3d& operator*=(const float rhs)
	{
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		return *this;
	}

	vec3d& operator/=(const float rhs)
	{
		this->x /= rhs;
		this->y /= rhs;
		this->z /= rhs;
		return *this;
	}

	vec3d operator+(const vec3d& rhs) const
	{
		return { this->x + rhs.x, this->y + rhs.y, this->z + rhs.z };
	}

	vec3d operator-(const vec3d& rhs) const
	{
		return { this->x - rhs.x, this->y - rhs.y, this->z - rhs.z };
	}

	vec3d operator*(const float rhs) const
	{
		return { this->x * rhs, this->y * rhs, this->z * rhs };
	}

	vec3d operator/(const float rhs) const
	{
		return { this->x / rhs, this->y / rhs, this->z / rhs };
	}

	vec3d operator*(const mat4x4f& rhs) const
	{
		vec3d r;
		r.x = this->x * rhs.mat[0][0] + this->y * rhs.mat[1][0] + this->z * rhs.mat[2][0] + this->w * rhs.mat[3][0];
		r.y = this->x * rhs.mat[0][1] + this->y * rhs.mat[1][1] + this->z * rhs.mat[2][1] + this->w * rhs.mat[3][1];
		r.z = this->x * rhs.mat[0][2] + this->y * rhs.mat[1][2] + this->z * rhs.mat[2][2] + this->w * rhs.mat[3][2];
		r.w = this->x * rhs.mat[0][3] + this->y * rhs.mat[1][3] + this->z * rhs.mat[2][3] + this->w * rhs.mat[3][3];
		return r;
	}

	#pragma endregion

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

		m_MeshCube.LoadFromObjectFile("res/obj/LowPolyTree.obj");

		// Projection Matrix
		m_MatProj = mat4x4f::CreateProjection(90, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.0f);

		return true;
	}

    bool OnUserUpdate(float fElapseTime) override
	{
		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

		//Set Rotation Matrices
		mat4x4f matRotZ, matRotX;
		m_fTheta += 1.0f * fElapseTime;

		matRotZ = mat4x4f::CreateRotationZ(m_fTheta * 0.5f);
		matRotX = mat4x4f::CreateRotationX(m_fTheta);

		//Set Translation Matrix
		mat4x4f matTrans;
		matTrans = mat4x4f::CreateTranslation(0.0f, -2.0f, 12.0f);

		//Set World Matrix
		mat4x4f matWorld;
		matWorld = mat4x4f::CreateIdentity();
		matWorld = matRotZ * matRotX;
		matWorld *= matTrans;

		std::vector<triangle> trianglesToRaster;

		//Draw Triangles
		for (const auto &tri : m_MeshCube.triangles)
		{
			//rotates the cube
			triangle triProjected, triTransformed;

			triTransformed.p[0] = tri.p[0] * matWorld;
			triTransformed.p[1] = tri.p[1] * matWorld;
			triTransformed.p[2] = tri.p[2] * matWorld;

			vec3d normal, vecA, vecB; //normal is the normal vector to the plane made with the 2 vectors a and b
			vecA = triTransformed.p[1] - triTransformed.p[0];
			vecB = triTransformed.p[2] - triTransformed.p[0];

			normal = vecA.CrossProduct(vecB); //cross prod between the 2 lines to get the normal
			normal = normal.Normalize(); //normalizing the normal

			//vector between the surface normal and the camera
			vec3d projLine = triTransformed.p[0] - m_vCamera;

			//draw triangle only if the normal is facing the cam so if its z is negative
			if (projLine.DotProduct(normal) < 0)
			{
				//Lighting
				vec3d directionalLight = { 0.0f, 0.0f, -1.0f }; //light facing the the cube
				directionalLight = directionalLight.Normalize();

				float dotProdNormLight = normal.DotProduct(directionalLight);

				CHAR_INFO c = GetColor(dotProdNormLight);
				triTransformed.col = c.Attributes;
				triTransformed.sym = c.Char.UnicodeChar;

				//projects the triangles from 3d to 2d
				triProjected.p[0] = triTransformed.p[0] * m_MatProj;
				triProjected.p[1] = triTransformed.p[1] * m_MatProj;
				triProjected.p[2] = triTransformed.p[2] * m_MatProj;

				// Scale into view, we normalize the coordinates into a 2d-kind-of space by dividing the vector by the depth z.
				triProjected.p[0] = triProjected.p[0] / triProjected.p[0].w;
				triProjected.p[1] = triProjected.p[1] / triProjected.p[1].w;
				triProjected.p[2] = triProjected.p[2] / triProjected.p[2].w;

				triProjected.col = triTransformed.col;
				triProjected.sym = triTransformed.sym;


				//Scale into view

				//will make range -1 / 1 to 0 / 2
				vec3d vOffsetView = { 1, 1, 0 };
				triProjected.p[0] += vOffsetView;
				triProjected.p[1] += vOffsetView;
				triProjected.p[2] += vOffsetView;

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
