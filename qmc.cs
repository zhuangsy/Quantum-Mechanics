/*********************************************************************
 * NOTE
 * COPYRIGHT Zhuang Siyuan 2014
 * This file is used for quantum-mechanics' calculation and expression
**********************************************************************/
using System;
using System.Text;
using System.Diagnostics;
using System.Numerics;

namespace qmc
{

    #region Complex Matrix
    [Serializable]
    public class mat : ICloneable
    {
        internal Complex[,] matrix;

        public enum VectorType { bra, ket }

        #region property
        public Complex this[int a, int b]
        {
            get { return matrix[a, b]; }
            set { matrix[a, b] = value; }
        }
        public int Height { get { return matrix.GetLength(0); } }
        public int Width { get { return matrix.GetLength(1); } }
        public bool Isbra { get { return Height == 1; } }
        public bool Isket { get { return Width == 1; } }
        public bool IsVector { get { return Isbra ^ Isket; } }
        public bool IsSquare { get { return Height == Width; } }
        public bool IsHermitate
        {
            get
            {
                if (!IsSquare) return false;
                for (int i = 0; i < Width; i++)
                    for (int j = i; j < Width; j++)
                        if (matrix[i, j] != Complex.Conjugate(matrix[j, i])) return false;
                return true;
            }
        }
        #endregion

        #region constructor
        public mat() { }
        public mat(Complex[,] m) { matrix = m; }
        public mat(int h, int w) { matrix = new Complex[h, w]; }
        public mat(int n) { matrix = new Complex[n, n]; }
        public static mat CreateBySize(mat m) { return new mat(m.Height, m.Width); }
        #endregion

        #region DefinedMatrixes
        public static mat Identity(int n)
        {
            var m = new mat(n);
            for (int i = 0; i < n; i++) m[i, i] = Complex.One;
            return m;
        }
        public static mat I
        {
            get { return new mat(new Complex[2, 2] { { Complex.One, Complex.Zero }, { Complex.Zero, Complex.One } }); }
        }
        public static mat PauliX
        {
            get
            {
                return new mat(new Complex[2, 2] { { 0, 1 }, { 1, 0 } });
            }
        }
        public static mat PauliY
        {
            get
            {
                return new mat(new Complex[2, 2] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } });
            }
        }
        public static mat PauliZ
        {
            get
            { return new mat(new Complex[2, 2] { { 1.0, 0.0 }, { 0.0, -1.0 } }); }
        }
        public static mat S
        {
            get { return new mat(new Complex[2, 2] { { 1.0, 0.0 }, { 0.0, Complex.ImaginaryOne } }); }
        }
        public static mat Hadamard
        {
            get
            {
                var t = new Complex(Math.Sqrt(0.5), 0);
                return new mat(new Complex[2, 2] { { t, t }, { t, -t } });
            }
        }
        public static mat CNOT
        {
            get
            {
                var m = new mat(4);
                m[0, 0] = Complex.One;
                m[1, 1] = Complex.One;
                m[2, 3] = Complex.One;
                m[3, 2] = Complex.One;
                return m;
            }
        }
        public static mat SWAP
        {
            get
            {
                var m = new mat(4);
                m[0, 0] = Complex.One;
                m[1, 2] = Complex.One;
                m[2, 1] = Complex.One;
                m[3, 3] = Complex.One;
                return m;
            }
        }
        #endregion

        #region operators
        public static bool CanMultiply(mat left, mat right) { return left.Width == right.Height; }
        public static bool IsSameInSize(mat left, mat right)
        { return (left.Width == right.Width) && (left.Height == right.Height); }
        public mat Transport()
        {
            var m = new mat(Width, Height);
            for (int i = 0; i < Height; i++)
                for (int j = 0; j < Width; j++)
                    m[j, i] = matrix[i, j];
            return m;
        }
        public static mat Transport(mat m) { return m.Transport(); }
        public mat Conjudgate()
        {
            var m = mat.CreateBySize(this);
            for (int i = 0; i < Height; i++)
                for (int j = 0; j < Width; j++)
                    m[i, j] = Complex.Conjugate(matrix[i, j]);
            return m;
        }
        public static mat Conjudgate(mat m) { return m.Conjudgate(); }
        public mat ConjudgateTransport()
        {
            var m = new mat(Width, Height);
            for (int i = 0; i < Height; i++)
                for (int j = 0; j < Width; j++)
                    m[j, i] = Complex.Conjugate(matrix[i, j]);
            return m;
        }
        public static mat ConjudgateTransport(mat m) { return m.ConjudgateTransport(); }

        #region Inner Product : 内积
        public static mat operator *(mat left, mat right)
        {
            if (!CanMultiply(left, right)) throw CannotMultiply;
            var m = new mat(left.Height, right.Width);

            for (int i = 0; i < left.Height; i++)
                for (int j = 0; j < right.Width; j++)
                    for (int k = 0; k < left.Width; k++)
                        m[i, j] += left[i, k] * right[k, j];
            return m;
        }
        public static ket operator *(mat left, ket right)
        {
            if (!CanMultiply(left, right)) throw CannotMultiply;
            var m = new mat(left.Height, 1);

            for (int i = 0; i < left.Height; i++)
                for (int j = 0; j < right.Width; j++)
                    for (int k = 0; k < left.Width; k++)
                        m[i, j] += left[i, k] * right[k, j];
            return m.Toket();
        }
        public static bra operator *(bra left, mat right)
        {
            if (!CanMultiply(left, right)) throw CannotMultiply;
            var m = new mat(1, right.Width);
            for (int i = 0; i < left.Height; i++)
                for (int j = 0; j < right.Width; j++)
                    for (int k = 0; k < left.Width; k++)
                        m[i, j] += left[i, k] * right[k, j];
            return m.Tobra();
        }
        #endregion

        /// <summary>plus 矩阵加法</summary>
        public static mat operator +(mat left, mat right)
        {
            if (!IsSameInSize(left, right)) throw SizeNotEqual;
            var m = CreateBySize(left);
            for (int i = 0; i < left.Height; i++)
                for (int j = 0; j < left.Width; j++)
                    m[i, j] = left[i, j] + right[i, j];
            return m;
        }

        #region Numberic Multiply：数乘
        public static mat operator *(mat left, Complex right)
        {
            var m = CreateBySize(left);
            for (int i = 0; i < left.Height; i++)
                for (int j = 0; j < left.Width; j++)
                    m[i, j] = left[i, j] * right;
            return m;
        }
        public static mat operator *(Complex left, mat right) { return right * left; }
        public static mat operator *(double left, mat right) { return right * left; }
        public static mat operator *(mat left, double right) { return ((Complex)right) * left; }
        #endregion

        #region Tensor product : 张量积
        /// <summary>tensor product 张量积</summary>
        public static mat operator ^(mat left, mat right)
        {
            var m = new mat(left.Height * right.Height, left.Width * right.Width);
            for (int ly = 0; ly < left.Height; ly++)
                for (int lx = 0; lx < left.Width; lx++)
                    for (int ry = 0; ry < right.Height; ry++)
                        for (int rx = 0; rx < right.Width; rx++)
                            m[ly * right.Height + ry, lx * right.Width + rx] = left[ly, lx] * right[ry, rx];
            return m;
        }
        public static mat operator <<(mat left, int n)
        {
            var m = left.Clone();
            for (int i = 0; i < n; i++) m = m ^ I;
            return m;
        }
        public static mat operator >>(mat left, int n)
        {
            var m = left.Clone();
            for (int i = 0; i < n; i++) m = I ^ m;
            return m;
        }
        #endregion

        #endregion

        #region TypeConvertion
        public virtual ket Toket()
        {
            if (Isket)
            {
                var k = new ket();
                k.matrix = (Complex[,])matrix.Clone();
                return k;
            }
            else
            {
                throw new InvalidCastException();
            }
        }
        public virtual bra Tobra()
        {
            if (Isbra)
            {
                var b = new bra();
                b.matrix = (Complex[,])matrix.Clone();
                return b;
            }
            else
            {
                throw new InvalidCastException();
            }
        }
        #endregion

        #region Utils
        public bool PossibilityCheck(Complex c)
        {
            return (c.Imaginary == 0) && (c.Real >= 0) && (c.Real <= 1);
        }
        static Random r = new Random();
        static Complex GetRandomComplex() { return new Complex(r.NextDouble(), r.NextDouble()); }
        public static ket GetRandomUnitKet()
        {
            var k = ket.Create(new Complex[] { GetRandomComplex(), GetRandomComplex() });
            k.Unitize();
            return k;
        }
        public static bra GetRandomUnitBra()
        {
            var b = bra.Create(new Complex[] { GetRandomComplex(), GetRandomComplex() });
            b.Unitize();
            return b;
        }
        #endregion

        #region Projection
        public static mat PZero
        {
            get { return new mat(new Complex[2, 2] { { 1.0, 0.0 }, { 0.0, 0.0 } }); }
        }
        public static mat POne
        {
            get { return new mat(new Complex[2, 2] { { 0.0, 0.0 }, { 0.0, 1.0 } }); }
        }
        /// <summary>
        /// Create_Single_Measuring_Operation_0
        /// </summary>
        public mat CSMO0(int index, int length)
        {
            var p = new mat(1);
            p[0, 0] = 1;
            for (int i = 0; i < length; i++)
                if (i != index) p ^= I; else p ^= PZero;
            return p;
        }
        public mat CSMO1(int index, int length)
        {
            var p = new mat(1);
            p[0, 0] = 1;
            for (int i = 0; i < length; i++)
                if (i != index) p ^= I; else p ^= POne;
            return p;
        }

        public double MeasureD(int index)
        {
            var bits = (int)Math.Floor(Math.Log(this.Width, 2));
            var r = AverageValueByDM(CSMO0(index, bits));
            if (!PossibilityCheck(r)) throw ErrProbability;
            return r.Real;
        }
        #endregion

        #region Density Matrix

        public Complex Trace()
        {
            if (!IsSquare) throw SquareMatrixNeeded;
            var tr = new Complex();
            for (int i = 0; i < Width; i++) tr += matrix[i, i];
            return tr;
        }
        public static Complex Trace(mat m)
        {
            return m.Trace();
        }
        public bool IsDensityMatrix//应该是必要条件
        {
            get { return IsHermitate && ((Trace() - Complex.One).Magnitude <= 0.000000001); }//考虑误差
        }

        #endregion

        #region AverageValue
        public Complex AverageValueByDM(mat m) //The average value(s) of an observable should be real
        {
            if (IsDensityMatrix) return (Trace(this * m)); else throw NotADensityMatrix;
        }

        #endregion

        #region expression

        internal string Dec2Bin(int n, int bits)
        {
            var s = string.Empty;
            var t = n;
            for (int i = 0; i < bits; i++)
            {
                s = (t % 2).ToString() + s;
                t = t / 2;
            }
            if (t > 0) Debug.Fail("Bits not enough.");
            return s;
        }

        public override string ToString()
        {
            var s = new StringBuilder();
            for (int i = 0; i < Height; i++)
            {
                for (int j = 0; j < Width; j++)
                {
                    s.Append("|" + matrix[i, j].ToString());
                }
                s.AppendLine("|");
            }
            return s.ToString();
        }
        #endregion

        #region Clone
        public mat Clone() { return new mat((Complex[,])matrix.Clone()); }
        object ICloneable.Clone() { return new mat((Complex[,])matrix.Clone()); }
        #endregion

        #region Error Processing
        internal protected static Exception NotAVector
        {
            get { return new Exception("Not a vector"); }
        }
        internal protected static Exception SquareMatrixNeeded
        {
            get { return new Exception("Square Matrix Needed"); }
        }
        internal protected static Exception NotADensityMatrix
        {
            get { return new Exception("Not a density matrix"); }
        }
        internal protected static Exception ErrProbability
        {
            get { return new Exception("Probability should be real,positive and less than one"); }
        }
        internal protected static Exception NotSpecificVectorType
        {
            get { return new Exception("VectorType should be either ket or bra"); }
        }
        internal protected static Exception SizeNotEqual
        {
            get { return new Exception("Size should be equal"); }
        }
        internal protected static Exception CannotMultiply
        {
            get { return new Exception("Not fit for multiply"); }
        }
        #endregion
    }
    #endregion

    #region Complex Vector
    //NOTE:|1000000=0〉
    public class ket : StateVector
    {
        public override Complex this[int index]
        {
            get { return matrix[index, 0]; }
            set { matrix[index, 0] = value; }
        }
        public override int VectorLength { get { return Height; } }
        public static ket Create(int n)
        {
            var k = new ket();
            k.matrix = new Complex[n, 1];
            return k;
        }
        public static ket Create(Complex[] v)
        {
            var k = new ket();
            k.matrix = new Complex[v.Length, 1];
            for (int i = 0; i < v.Length; i++) k.matrix[i, 0] = v[i];
            return k;
        }

        public override mat CreateDensityMatrix() { return this ^ ConjudgateTransport(); }
        public static mat CreateDensityMatrix(ket k) { return k.CreateDensityMatrix(); }
        public static Complex operator *(bra left, ket right)
        {
            if (!CanMultiply(left, right)) throw CannotMultiply;
            var c = Complex.Zero;
            for (int i = 0; i < left.VectorLength; i++) c += left[i] * right[i];
            return c;
        }
        public static ket operator ^(ket left, ket right) { return (((mat)left) ^ ((mat)right)).Toket(); }
        #region Numberic Multiply
        public static ket operator *(Complex left, ket right)
        {
            var k = Create(right.VectorLength);
            for (int i = 0; i < right.VectorLength; i++) { k[i] = right[i] * left; }
            return k;
        }
        public static ket operator *(ket left, Complex right) { return right * left; }
        public static ket operator *(double left, ket right) { return (Complex)left * right; }
        public static ket operator *(ket left, double right) { return right * left; }
        #endregion

        public override bra Tobra() { return ConjudgateTransport().Tobra(); }
        public static ket ket0 { get { return ket.Create(zero); } }
        public static ket ket1 { get { return ket.Create(one); } }

        public override Complex AverageValueUnderCertainState(mat m) { return ((Tobra() * m) * this); }
        public override string Vectorprint
        {
            get
            {
                var bits = VectorBits;
                var sb = new StringBuilder();
                for (int i = 0; i < VectorLength; i++)
                    sb.AppendLine(this[i].ToString() + "|" + Dec2Bin(i, bits) + ">");
                return sb.ToString();
            }
        }
    }

    public class bra : StateVector
    {
        public override Complex this[int index]
        {
            get { return matrix[0, index]; }
            set { matrix[0, index] = value; }
        }
        public override int VectorLength { get { return Width; } }
        public static bra Create(Complex[] v)
        {
            var b = new bra();
            b.matrix = new Complex[1, v.Length];
            for (int i = 0; i < v.Length; i++) b.matrix[0, i] = v[i];
            return b;
        }
        public static bra Create(int n)
        {
            var b = new bra();
            b.matrix = new Complex[n, 1];
            return b;
        }
        public override mat CreateDensityMatrix() { return ConjudgateTransport() ^ this; }
        public static mat CreateDensityMatrix(bra b) { return b.CreateDensityMatrix(); }

        #region Numberic Multiply
        public static bra operator *(Complex left, bra right)
        {
            var b = Create(right.VectorLength);
            for (int i = 0; i < right.VectorLength; i++) { b[i] = right[i] * left; }
            return b;
        }
        public static bra operator *(bra left, Complex right) { return right * left; }
        public static bra operator *(double left, bra right) { return (Complex)left * right; }
        public static bra operator *(bra left, double right) { return right * left; }
        #endregion
        public static bra operator ^(bra left, bra right) { return (((mat)left) ^ ((mat)right)).Tobra(); }

        public static bra bra0 { get { return bra.Create(zero); } }
        public static bra bra1 { get { return bra.Create(one); } }

        public override ket Toket() { return ConjudgateTransport().Toket(); }
        public override Complex AverageValueUnderCertainState(mat m) { return (this * m * Toket()); }

        public override string Vectorprint
        {
            get
            {
                var bits = VectorBits;
                var sb = new StringBuilder();
                for (int i = 0; i < VectorLength; i++)
                    sb.AppendLine(this[i].ToString() + "<" + Dec2Bin(i, bits) + "|");
                return sb.ToString();
            }
        }
    }

    public abstract class StateVector : mat
    {
        public static Complex[] zero { get { return new Complex[] { Complex.One, Complex.Zero }; } }
        public static Complex[] one { get { return new Complex[] { Complex.Zero, Complex.One }; } }
        public abstract Complex this[int index] { get; set; }
        public abstract int VectorLength { get; }
        public int VectorBits
        {
            get { return (int)Math.Floor(Math.Log((double)(VectorLength - 1), 2)) + 1; }
        }
        public void Unitize()
        {
            var c = 0.0;
            for (int i = 0; i < VectorLength; i++)
                c += this[i].Real * this[i].Real + this[i].Imaginary * this[i].Imaginary;
            c = Math.Sqrt(c);
            for (int i = 0; i < VectorLength; i++) this[i] = this[i] / c;
        }
        //The possiblity of measuring zero
        public double MeasureV(int index)
        {
            var bits = (int)Math.Floor(Math.Log(VectorLength, 2));
            var r = AverageValueUnderCertainState(CSMO0(index, bits));
            if (!PossibilityCheck(r)) throw ErrProbability;
            return r.Real;
        }
        public abstract mat CreateDensityMatrix();
        //The average value(s) of an observable should be real
        public abstract Complex AverageValueUnderCertainState(mat m);
        public abstract string Vectorprint { get; }
    }

    #endregion

}
