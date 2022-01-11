/*
 * SPDX-FileCopyrightText: 2020 Kevin Höllring for PULS Group <kevin.hoellring@fau.de>
 *
 * SPDX-License-Identifier: MIT
 * 
 * Copyright (c) 2020 Kevin Höllring for PULS Group <kevin.hoellring@fau.de>
 * 
 * Authors: 2020 Kevin Höllring <kevin.hoellring@fau.de>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy 
 * of this software and associated documentation files (the “Software”), to deal 
 * in the Software without restriction, including without limitation the rights 
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
 * copies of the Software, and to permit persons to whom the Software is 
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice (including the next paragraph) 
 * shall be included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE. 
 */

#ifndef __MOLECULAR_GRAPH_VEC_H__
#define __MOLECULAR_GRAPH_VEC_H__

/**
 * @brief Simple geometric 3d vector class
 * 
 * @tparam T 
 */
template <typename T>
class vec
{
public:
    union
    {
        struct
        {
            T x, y, z;
        };
        T pos[3];
    };

    vec() : x(T()), y(T()), z(T()) {}
    vec(const vec<T> &other) : x(other.x), y(other.y), z(other.z) {}

    vec<T> &operator=(const vec<T> &other)
    {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }

    vec<T> operator+(const vec<T> &other) const
    {
        vec<T> res(*this);
        res += other;
        return res;
    }

    vec<T> operator-(const vec<T> &other) const
    {
        vec<T> res(*this);
        res -= other;
        return res;
    }

    T operator*(const vec<T> &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

    vec<T> operator*(T scalar) const
    {
        vec<T> res(*this);
        res *= scalar;
        return res;
    }

    vec<T> operator/(T scalar) const
    {
        vec<T> res(*this);
        res /= scalar;
        return res;
    }

    vec<T> &operator+=(const vec<T> &other)
    {
        for (size_t i = 0; i < 3; i++)
        {
            pos[i] += other.pos[i];
        }
        return *this;
    }

    vec<T> &operator-=(const vec<T> &other)
    {
        for (size_t i = 0; i < 3; i++)
        {
            pos[i] -= other.pos[i];
        }
        return *this;
    }

    vec<T> &operator*=(T scalar)
    {
        for (size_t i = 0; i < 3; i++)
        {
            pos[i] *= scalar;
        }
        return *this;
    }

    vec<T> &operator/=(T scalar)
    {
        for (size_t i = 0; i < 3; i++)
        {
            pos[i] /= scalar;
        }
        return *this;
    }

    T operator[](size_t i) const
    {
        return pos[i];
    }

    T &operator[](size_t i)
    {
        return pos[i];
    }

    T norm() const
    {
        return std::sqrt(norm2());
    }

    T norm2() const
    {
        T res = T();
        for (size_t i = 0; i < 3; i++)
        {
            res += pos[i] * pos[i];
        }
        return res;
    }
};

/**
 * @brief Operator to allow for both orders of scalar multiplication
 * 
 * @tparam T 
 * @param scalar 
 * @param other 
 * @return vec<T> 
 */
template <typename T>
vec<T> operator*(T scalar, const vec<T> &other)
{
    vec<T> res(other);
    res *= scalar;
    return res;
}

/**
 * @brief Retrieve the basis coefficients of vector @p pos with respect to a triclinic @p basis 
 * 
 * @tparam T 
 * @param pos 
 * @param basis 
 * @return vec<T> 
 */
template <typename T>
vec<T> decompose(const vec<T> &pos, const std::vector<vec<T>> &basis)
{
    vec<T> cpy(pos);
    vec<T> res;

    for (int i = 2; i >= 0; i--)
    {
        res[i] = cpy[i] / basis[i][i];
        cpy -= basis[i] * res[i];
    }
    return res;
}

/**
 * @brief Removes pbc from basis coefficients to limit them to the interval [-0.5, 0.5] respectively
 * 
 * @tparam T 
 * @param coeff 
 * @return vec<T> 
 */
template <typename T>
vec<T> normalize_basis_coefficients(const vec<T> &coeff)
{
    vec<T> res(coeff);
    for (size_t c = 0; c < 3; c++)
    {
        // This step removes pbc for coefficients being between [0,1]
        res[c] -= std::floor(res[c]);
        // This step takes the coefficients to [-0.5, 0.5]
        if (res[c] >= 0.5)
        {
            res[c] -= 1.0;
        }
    }
    return res;
}
#endif