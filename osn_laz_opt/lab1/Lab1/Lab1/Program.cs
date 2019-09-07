using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lab1
{
    class Program
    {
        public static matrica umnozh(matrica m1_, matrica m2_)
        {
            double pa = m1_.a * m2_.a + m1_.b * m2_.c;
            double pb = m1_.a * m2_.b + m1_.b * m2_.d;
            double pc = m1_.c * m2_.a + m1_.d * m2_.c;
            double pd = m1_.c * m2_.b + m1_.d * m2_.d;
            matrica p_ = new matrica(pa, pb, pc, pd);
            return p_;
        }
        public static double radius(matrica XX, double YZR, double YZ0, double YLAMBDA)
        {
            double o1 = -1 * YZR * XX.a;
            double o2 = XX.a * YZ0 + XX.b;
            double o3 = -1 * YZR * XX.c;
            double o4 = XX.c * YZ0 + XX.d;
            double qRe = ((o2 * o4) + (o1 * o3)) / ((o4 * o4) + (o3 * o3));
            double qIm = ((o1 * o4) - (o2 * o3)) / ((o4 * o4) + (o3 * o3));
            double radius_ = Math.Sqrt(-1 * (YLAMBDA / Math.PI) * qIm * (1 + ((qRe / qIm) * (qRe / qIm))));
            return radius_;
        }
        static void Main(string[] args)
        {
            // variant 6
            double lambda = 0.002570;
            double R1 = 400;
            double R2 = 400;
            double f = 200;
            double L = 500;
            double L1 = 200;
            double Dfiber = 0.6;
            double NAfiber = 0.22;

            double L2 = L - L1;
            double L_ = L1 + L2 - ((L1 * L2) / f);
            double g1 = 1 - (L2 / f) - (L_ / R1);
            double g2 = 1 - (L1 / f) - (L_ / R2);

            double w1 = Math.Sqrt((lambda * L_ / Math.PI) * Math.Sqrt(g2 / (g1 * (1 - (g1 * g2)))));
            double w2 = Math.Sqrt((lambda * L_ / Math.PI) * Math.Sqrt(g1 / (g2 * (1 - (g1 * g2)))));

            double z0 = R2 / (1 + ((R2 / (Math.PI * w2 * w2 / lambda)) * (R2 / (Math.PI * w2 * w2 / lambda))));
            double w0 = Math.Sqrt(((w2 * w2) - Math.Sqrt((w2 * w2 * w2 * w2) - (4 * ((lambda * z0 / Math.PI) * (lambda * z0 / Math.PI))))) / 2);

            double thetamax = Math.Asin(NAfiber);

            double zR = Math.PI * w0 * w0 / lambda;

            double w = 0;
            double theta = 0;
            double dw = 10;
            double dtheta = 10;
            double mw = 0;
            double mtheta = 0;
            double mf1 = 0;
            double mLObj = 0;
            double mf2 = 0;

            matrica st1 = new matrica(0, 0, 0, 0);
            matrica st2 = new matrica(0, 0, 0, 0);
            matrica st3 = new matrica(0, 0, 0, 0);
            matrica st4 = new matrica(0, 0, 0, 0);
            matrica st5 = new matrica(0, 0, 0, 0);
            matrica finalM = new matrica(0, 0, 0, 0);

            for (double f1 = 200; f1 <= 1200; f1 += 10)
                for (double LObj = 10; LObj <= 1010; LObj += 10)
                    for (double f2 = 10; f2 <= 1010; f2 += 10)
                    {
                        st1.a = 1; st1.b = f1; st1.c = 0; st1.d = 1;
                        st2.a = 1; st2.b = 0; st2.c = (-1 / f1); st2.d = 1;
                        st3.a = 1; st3.b = LObj; st3.c = 0; st3.d = 1;
                        st4.a = 1; st4.b = 0; st4.c = (-1 / f2); st4.d = 1;
                        st5.a = 1; st5.b = f2; st5.c = 0; st5.d = 1;
                        finalM = umnozh((umnozh((umnozh((umnozh(st5, st4)), st3)), st2)), st1);
                        w = radius(finalM, zR, z0, lambda);
                        theta = (1.22 * lambda) / (2 * w);
                        if ((w <= (Dfiber / 2)) & (theta <= thetamax) & ((Math.Abs((Dfiber / 2) - w)) < dw))
                        {
                            dw = Math.Abs((Dfiber / 2) - w);
                            dtheta = Math.Abs(thetamax - theta);
                            mw = w;
                            mtheta = theta;
                            mf1 = f1;
                            mLObj = LObj;
                            mf2 = f2;
                        }
                    }

            double mmf1 = mf1;
            double mmLObj = mLObj;
            double mmf2 = mf2;

            for (double f1 = mmf1 - 10; f1 <= mmf1 + 10; f1 += 0.5)
                for (double LObj = mmLObj - 10; LObj <= mmLObj + 10; LObj += 0.5)
                    for (double f2 = mmf2 - 10; f2 <= mmf2 + 10; f2 += 0.5)
                    {
                        st1.a = 1; st1.b = f1; st1.c = 0; st1.d = 1;
                        st2.a = 1; st2.b = 0; st2.c = (-1 / f1); st2.d = 1;
                        st3.a = 1; st3.b = LObj; st3.c = 0; st3.d = 1;
                        st4.a = 1; st4.b = 0; st4.c = (-1 / f2); st4.d = 1;
                        st5.a = 1; st5.b = f2; st5.c = 0; st5.d = 1;
                        finalM = umnozh((umnozh((umnozh((umnozh(st5, st4)), st3)), st2)), st1);
                        w = radius(finalM, zR, z0, lambda);
                        theta = (1.22 * lambda) / (2 * w);
                        if ((w <= (Dfiber / 2)) & (theta <= thetamax) & ((Math.Abs((Dfiber / 2) - w)) < dw))
                        {
                            dw = Math.Abs((Dfiber / 2) - w);
                            dtheta = Math.Abs(thetamax - theta);
                            mw = w;
                            mtheta = theta;
                            mf1 = f1;
                            mLObj = LObj;
                            mf2 = f2;
                        }
                    }

            int zN = (int)mf1 + (int)mLObj + (int)mf2 + 1;
            double[] wz = new double[zN];
            for (int i = 0; i <= (int)mf1; i++)
            {
                st1.a = 1; st1.b = i; st1.c = 0; st1.d = 1;
                finalM = st1;
                wz[i] = radius(finalM, zR, z0, lambda);
            }
            for (int i = (int)mf1 + 1; i <= (int)mf1 + (int)mLObj; i++)
            {
                st1.a = 1; st1.b = mf1; st1.c = 0; st1.d = 1;
                st2.a = 1; st2.b = 0; st2.c = (-1 / mf1); st2.d = 1;
                st3.a = 1; st3.b = i - mf1; st3.c = 0; st3.d = 1;
                finalM = umnozh((umnozh(st3, st2)), st1);
                wz[i] = radius(finalM, zR, z0, lambda);
            }
            for (int i = (int)mf1 + (int)mLObj + 1; i <= zN - 1; i++)
            {
                st1.a = 1; st1.b = mf1; st1.c = 0; st1.d = 1;
                st2.a = 1; st2.b = 0; st2.c = (-1 / mf1); st2.d = 1;
                st3.a = 1; st3.b = mLObj; st3.c = 0; st3.d = 1;
                st4.a = 1; st4.b = 0; st4.c = (-1 / mf2); st4.d = 1;
                st5.a = 1; st5.b = i - mf1 - mLObj; st5.c = 0; st5.d = 1;
                finalM = umnozh((umnozh((umnozh((umnozh(st5, st4)), st3)), st2)), st1);
                wz[i] = radius(finalM, zR, z0, lambda);
            }




            Console.WriteLine("g1 = " + g1.ToString());
            Console.WriteLine("g2 = " + g2.ToString());
            Console.WriteLine();
            Console.WriteLine("w0 = " + w0.ToString());
            Console.WriteLine("z0 = " + z0.ToString());
            Console.WriteLine();
            Console.WriteLine();
            Console.WriteLine("Dfiber = " + Dfiber.ToString());
            Console.WriteLine("2w = " + (2 * mw).ToString());
            Console.WriteLine();
            Console.WriteLine("thetamax = " + thetamax.ToString());
            Console.WriteLine("theta = " + mtheta.ToString());
            Console.WriteLine();
            Console.WriteLine();
            Console.WriteLine("f1 = " + mf1.ToString());
            Console.WriteLine();
            Console.WriteLine("LObj = " + mLObj.ToString());
            Console.WriteLine();
            Console.WriteLine("f2 = " + mf2.ToString());

            Console.WriteLine();
            Console.WriteLine();
            Console.WriteLine();
            for (int i = (int)z0 + 1; i <= zN - 1; i++)
            {
                Console.WriteLine((int)(i - z0));
            }
            Console.WriteLine();
            for (int i = (int)z0 + 1; i <= zN - 1; i++)
            {
                Console.WriteLine(wz[i]);
            }


            Console.ReadLine();
        }
    }
    public class matrica
    {
        private double a_;
        public double a
        {
            get { return a_; }
            set { a_ = value; }
        }
        private double b_;
        public double b
        {
            get { return b_; }
            set { b_ = value; }
        }
        private double c_;
        public double c
        {
            get { return c_; }
            set { c_ = value; }
        }
        private double d_;
        public double d
        {
            get { return d_; }
            set { d_ = value; }
        }
        public matrica(double AA, double BB, double CC, double DD)
        {
            a = AA;
            b = BB;
            c = CC;
            d = DD;
        }
    }
}
