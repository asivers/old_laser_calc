using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApplication1
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
        static void Main(string[] args)
        {
            //задаем переменные
            double qv = 8;
            double Tc = 293;
            double k = 0.14;
            double a = 0.5;
            double Rae = 1;
            double Lae = 8;
            double P = 0.0000082;
            double Q = 0.0000017;
            double n0 = 1.8197;
            double R1 = -200;
            double R2 = -200;
            double Lc = 12;
            double Lc1 = 6;

            //расчет температуры
            double To = Tc + ((qv * Rae * Rae) / (4 * k)) + ((qv * Rae) / (2 * a));
            double dT = (qv * Rae * Rae) / (4 * k);

            //цикл для построения графика зависимости температуры от радиуса
            double[] T = new double[(int)(100 * Rae) + 1];
            for (int i = 0; i <= (int)(100 * Rae); i++)
            {
                T[i] = To - (dT * Math.Pow(Rae, -2) * i * i / 10000);
            }

            //расчет фокусного расстояния
            double fr = (Rae * Rae) / (2 * (P + (Q / 2)) * dT * Lae);
            double ff = (Rae * Rae) / (2 * (P - (Q / 2)) * dT * Lae);
            double f = (fr + ff) / 2;

            //задаем матрицы
            matrica tonkx = new matrica(0, 0, 0, 0);
            matrica prex = new matrica(0, 0, 0, 0);
            matrica aftx = new matrica(0, 0, 0, 0);
            matrica x = new matrica(0, 0, 0, 0);

            //считаем матрицу для пути в одну сторону для тонкой линзы
            tonkx.a = 1; tonkx.b = 0; tonkx.c = (-1 / f); tonkx.d = 1;
            prex.a = 1; prex.b = 0.5 * Lae / n0; prex.c = 0; prex.d = 1;
            aftx.a = 1; aftx.b = 0.5 * Lae / n0; aftx.c = 0; aftx.d = 1;
            x = umnozh((umnozh(aftx, tonkx)), prex);

            //считаем npov на поверхности и n2 (не показатель преломления, размерность см^-1)
            double npov = n0 + (P * dT / 2) - ((P + (Q / 2)) * dT * 1);
            double n2 = 2 * (n0 - npov) / (Rae * Rae);

            //упрощаем расчет для "толстой" линзы
            matrica y = new matrica(0, 0, 0, 0);
            double skob = Lae * Math.Sqrt(n2 / n0);
            double kor = Math.Sqrt(n2 * n0);

            //считаем матрицу для пути в одну сторону для "толстой" линзы
            y.a = Math.Cos(skob);
            y.b = (1 / kor) * Math.Sin(skob);
            y.c = (-1 * kor) * Math.Sin(skob);
            y.d = Math.Cos(skob);

            //задаем матрицы для полного обхода резонатора
            double Lc2 = Lc - Lc1;
            matrica z1 = new matrica(0, 0, 0, 0);
            matrica z2 = new matrica(0, 0, 0, 0);
            matrica pre = new matrica(0, 0, 0, 0);
            matrica aft = new matrica(0, 0, 0, 0);
            matrica polntonk = new matrica(0, 0, 0, 0);
            z1.a = 1; z1.b = 0; z1.c = (-2 / R1); z1.d = 1;
            z2.a = 1; z2.b = 0; z2.c = (-2 / R2); z2.d = 1;
            pre.a = 1; pre.b = Lc1 - (0.5 * Lae); pre.c = 0; pre.d = 1;
            aft.a = 1; aft.b = Lc2 - (0.5 * Lae); aft.c = 0; aft.d = 1;

            //считаем матрицу для полного обхода резонатора
            polntonk = umnozh(z1, pre);
            polntonk = umnozh(polntonk, prex);
            polntonk = umnozh(polntonk, tonkx);
            polntonk = umnozh(polntonk, aftx);
            polntonk = umnozh(polntonk, aft);
            polntonk = umnozh(polntonk, z2);
            polntonk = umnozh(polntonk, aft);
            polntonk = umnozh(polntonk, aftx);
            polntonk = umnozh(polntonk, tonkx);
            polntonk = umnozh(polntonk, prex);
            polntonk = umnozh(polntonk, pre);

            //проверяем неустойчивость резонатора (>1 неустойчивый)
            double prov = (polntonk.a + polntonk.d) / 2;

            //задаем матрицы для полного обхода эквивалентного резонатора
            matrica plz1 = new matrica(0, 0, 0, 0);
            matrica tl1 = new matrica(0, 0, 0, 0);
            matrica plz2 = new matrica(0, 0, 0, 0);
            matrica tl2 = new matrica(0, 0, 0, 0);
            matrica polnekv = new matrica(0, 0, 0, 0);
            plz1.a = 1; plz1.b = 0; plz1.c = 0; plz1.d = 1;
            tl1.a = 1; tl1.b = 0; tl1.c = (-1 / R1); tl1.d = 1;
            plz2.a = 1; plz2.b = 0; plz2.c = 0; plz2.d = 1;
            tl2.a = 1; tl2.b = 0; tl2.c = (-1 / R2); tl2.d = 1;

            //считаем матрицу для полного обхода эквивалентного резонатора
            polnekv = umnozh(tl1, plz1);
            polnekv = umnozh(polnekv, tl1);
            polnekv = umnozh(polnekv, pre);
            polnekv = umnozh(polnekv, prex);
            polnekv = umnozh(polnekv, tonkx);
            polnekv = umnozh(polnekv, aftx);
            polnekv = umnozh(polnekv, aft);
            polnekv = umnozh(polnekv, tl2);
            polnekv = umnozh(polnekv, plz2);
            polnekv = umnozh(polnekv, tl2);
            polnekv = umnozh(polnekv, aft);
            polnekv = umnozh(polnekv, aftx);
            polnekv = umnozh(polnekv, tonkx);
            polnekv = umnozh(polnekv, prex);
            polnekv = umnozh(polnekv, pre);

            //проверяем неустойчивость эквивалентного резонатора (>1 неустойчивый)
            double provekv = polnekv.a * polnekv.d;

            ////проверяем неустойчивость через g1, g2 (>1 неустойчивый)
            double L_ = Lc1 + Lc2 - (Lc1 * Lc2 / f);
            double g1 = 1 - (Lc2 / f) - (L_ / R1);
            double g2 = 1 - (Lc1 / f) - (L_ / R2);
            double provg = g1 * g2;

            //определяем R и mu для усредненного f
            double Rtrue = Math.Pow(R2, -1) + Math.Sqrt((polnekv.c * polnekv.d) / (polnekv.a * polnekv.b));
            Rtrue = Math.Pow(Rtrue, -1);
            double Rfalse = Math.Pow(R2, -1) - Math.Sqrt((polnekv.c * polnekv.d) / (polnekv.a * polnekv.b));
            Rfalse = Math.Pow(Rfalse, -1);
            double mutrue = (polnekv.a * polnekv.d) + (polnekv.b * polnekv.c) + 2 * Math.Sqrt(polnekv.a * polnekv.b * polnekv.c * polnekv.d);
            double mufalse = (polnekv.a * polnekv.d) + (polnekv.b * polnekv.c) - 2 * Math.Sqrt(polnekv.a * polnekv.b * polnekv.c * polnekv.d);

            //задаем матрицы для вычисления матрицы обхода резонатора для fr и ff
            matrica tonkxr = new matrica(0, 0, 0, 0);
            matrica tonkxf = new matrica(0, 0, 0, 0);
            tonkxr.a = 1; tonkxr.b = 0; tonkxr.c = (-1 / fr); tonkxr.d = 1;
            tonkxf.a = 1; tonkxf.b = 0; tonkxf.c = (-1 / ff); tonkxf.d = 1;
            matrica polnekvr = new matrica(0, 0, 0, 0);
            matrica polnekvf = new matrica(0, 0, 0, 0);

            //считаем R и mu для fr
            polnekvr = umnozh(tl1, plz1);
            polnekvr = umnozh(polnekvr, tl1);
            polnekvr = umnozh(polnekvr, pre);
            polnekvr = umnozh(polnekvr, prex);
            polnekvr = umnozh(polnekvr, tonkxr);
            polnekvr = umnozh(polnekvr, aftx);
            polnekvr = umnozh(polnekvr, aft);
            polnekvr = umnozh(polnekvr, tl2);
            polnekvr = umnozh(polnekvr, plz2);
            polnekvr = umnozh(polnekvr, tl2);
            polnekvr = umnozh(polnekvr, aft);
            polnekvr = umnozh(polnekvr, aftx);
            polnekvr = umnozh(polnekvr, tonkxr);
            polnekvr = umnozh(polnekvr, prex);
            polnekvr = umnozh(polnekvr, pre);
            double Rr = Math.Pow(R2, -1) + Math.Sqrt((polnekvr.c * polnekvr.d) / (polnekvr.a * polnekvr.b));
            Rr = Math.Pow(Rr, -1);
            double mur = (polnekvr.a * polnekvr.d) + (polnekvr.b * polnekvr.c) + 2 * Math.Sqrt(polnekvr.a * polnekvr.b * polnekvr.c * polnekvr.d);

            //считаем R и mu для ff
            polnekvf = umnozh(tl1, plz1);
            polnekvf = umnozh(polnekvf, tl1);
            polnekvf = umnozh(polnekvf, pre);
            polnekvf = umnozh(polnekvf, prex);
            polnekvf = umnozh(polnekvf, tonkxf);
            polnekvf = umnozh(polnekvf, aftx);
            polnekvf = umnozh(polnekvf, aft);
            polnekvf = umnozh(polnekvf, tl2);
            polnekvf = umnozh(polnekvf, plz2);
            polnekvf = umnozh(polnekvf, tl2);
            polnekvf = umnozh(polnekvf, aft);
            polnekvf = umnozh(polnekvf, aftx);
            polnekvf = umnozh(polnekvf, tonkxf);
            polnekvf = umnozh(polnekvf, prex);
            polnekvf = umnozh(polnekvf, pre);
            double Rf = Math.Pow(R2, -1) + Math.Sqrt((polnekvf.c * polnekvf.d) / (polnekvf.a * polnekvf.b));
            Rf = Math.Pow(Rf, -1);
            double muf = (polnekvf.a * polnekvf.d) + (polnekvf.b * polnekvf.c) + 2 * Math.Sqrt(polnekvf.a * polnekvf.b * polnekvf.c * polnekvf.d);

            //задаем переменные для подбора f тепловой линзы резонатора на границе устойчивости
            double fx = 0;
            double L_x = 0;
            double g1x = 0;
            double g2x = 0;
            double provgx = 0;
            double dprovgx = 100;
            double mfx = 0;
            double mg1x = 0;
            double mg2x = 0;

            //подбираем f тепловой линзы резонатора на границе устойчивости
            for (int i = 10; i <= 1000; i++)
            {
                fx = i;
                L_x = Lc1 + Lc2 - (Lc1 * Lc2 / fx);
                g1x = 1 - (Lc2 / fx) - (L_x / R1);
                g2x = 1 - (Lc1 / fx) - (L_x / R2);
                provgx = g1x * g2x;
                if ((Math.Abs(1 - provgx)) < dprovgx)
                {
                    dprovgx = Math.Abs(1 - provgx);
                    mfx = fx;
                    mg1x = g1x;
                    mg2x = g2x;
                }
            }

            Console.WriteLine("fr и ff тепловой линзы");
            Console.WriteLine(fr);
            Console.WriteLine(ff);
            Console.WriteLine();
            Console.WriteLine("f тепловой линзы");
            Console.WriteLine(f);
            Console.WriteLine();
            Console.WriteLine("тонкая линза: матрица пути в одну сторону");
            Console.WriteLine(x.a);
            Console.WriteLine(x.b);
            Console.WriteLine(x.c);
            Console.WriteLine(x.d);
            Console.WriteLine();
            Console.WriteLine("толстая линза: матрица пути в одну сторону");
            Console.WriteLine(y.a);
            Console.WriteLine(y.b);
            Console.WriteLine(y.c);
            Console.WriteLine(y.d);
            Console.WriteLine();
            Console.WriteLine("стандартный резонатор: матрица полного обхода");
            Console.WriteLine(polntonk.a);
            Console.WriteLine(polntonk.b);
            Console.WriteLine(polntonk.c);
            Console.WriteLine(polntonk.d);
            Console.WriteLine();
            Console.WriteLine("стандартный резонатор: критерий устойчивости");
            Console.WriteLine(prov);
            Console.WriteLine();
            Console.WriteLine("эквивалентный резонатор: матрица полного обхода");
            Console.WriteLine(polnekv.a);
            Console.WriteLine(polnekv.b);
            Console.WriteLine(polnekv.c);
            Console.WriteLine(polnekv.d);
            Console.WriteLine();
            Console.WriteLine("эквивалентный резонатор: критерий устойчивости");
            Console.WriteLine(provekv);
            Console.WriteLine();
            Console.WriteLine("критерий устойчивости через g1, g2");
            Console.WriteLine(provg);
            Console.WriteLine();
            Console.WriteLine("R и mu для усредненного f");
            Console.WriteLine(Rtrue);
            Console.WriteLine(mutrue);
            Console.WriteLine();
            Console.WriteLine("R и mu для fr");
            Console.WriteLine(Rr);
            Console.WriteLine(mur);
            Console.WriteLine();
            Console.WriteLine("R и mu для ff");
            Console.WriteLine(Rf);
            Console.WriteLine(muf);
            Console.WriteLine();
            Console.WriteLine("f, при котором резонатор находится на границе устойчивости");
            Console.WriteLine(mfx);
            Console.WriteLine();
            Console.WriteLine("g1 и g2");
            Console.WriteLine(mg1x);
            Console.WriteLine(mg2x);
            Console.WriteLine();
            Console.WriteLine();
            for (int i = 0; i <= (int)(100 * Rae); i++)
            {
                Console.WriteLine(i * 0.01);
            }
            Console.WriteLine();
            for (int i = 0; i <= (int)(100 * Rae); i++)
            {
                Console.WriteLine(T[i]);
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
