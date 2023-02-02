using DelaunatorSharp;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace TrianglifyCS
{
	/// <summary>
	/// Interaction logic for MainWindow.xaml
	/// </summary>
	public partial class MainWindow : Window
	{
		Random rand = new Random();

		public MainWindow()
		{
			InitializeComponent();
			Run();
		}

		public void Run()
		{
			var width = 1920;
			var height = 1080;
			var minX = 0;
			var maxX = minX + width;
			var minY = 0;
			var maxY = minY + height;
			var divX = 16;
			var divY = 9;
			var repeatX = 3;
			var repeatY = 3;
			var delta = 50;
			var deltaBright = 0.1;
			var baseBright = 0.92;

			var lineStart = new Complex(0, 0);
			var lineEnd = new Complex(1920, 1080);

			var colors = new List<(Color c, double mile)>();
			colors.Add((Color.FromRgb(247, 77, 35), 0.0));
			colors.Add((Color.FromRgb(240, 218, 26), 0.3));
			colors.Add((Color.FromRgb(70, 224, 43), 0.5));
			colors.Add((Color.FromRgb(240, 218, 26), 0.7));
			colors.Add((Color.FromRgb(240, 105, 10), 1.0));
			var line = new Geometry2d.Line2d(lineStart, lineEnd);

			this.Height = maxX - minX;
			this.Width = maxY - minY;

			var canvas = new Canvas();
			canvas.Width = repeatX * width;
			canvas.Height = repeatY * height;

			var points = GeneratePoints(minX, maxX, minY, maxY, divX, divY, repeatX, repeatY, delta);
			var delaunay = new Delaunator(points.Select(p => (IPoint)p).ToArray());
			var tri = delaunay.GetTriangles();

			var dic = new List<(double t, Color c)>();
			var polys = new List<Polygon>();
			foreach (var triangle in tri)
			{
				var poly = new Polygon();
				var g = new System.Windows.Point(triangle.Points.Select(p => p.X).Sum() / 3.0, triangle.Points.Select(p => p.Y).Sum() / 3.0);
				var pos = CalcPosition(g, line);
				var color = new Color();
				color = GetColor((pos + 1) % 1, colors, deltaBright, baseBright);
				dic.Add((pos, color));
				poly.Fill = new SolidColorBrush(color);
				foreach (var p in triangle.Points)
				{
					poly.Points.Add(new System.Windows.Point(p.X, p.Y));
				}
				canvas.Children.Add(poly);
				polys.Add(poly);
			}

			canvas.toImage($"image_{DateTime.Now.ToString("yyyyMMdd")}.png");
			this.Content = canvas;
		}

		public double Dist(System.Windows.Point p1, System.Windows.Point p2) => Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y));

		public Color GetColor(double t, List<(Color color, double mile)> colors, double delta, double baseBright)
		{
			for (int i = 0; i < colors.Count - 1; i++)
			{
				if (t < colors[i].mile || colors[i + 1].mile <= t) continue;
				var t1 = (t - colors[i].mile) / (colors[i + 1].mile - colors[i].mile);
				var c0 = colors[i].color;
				var c1 = colors[i + 1].color;
				var s = (rand.NextDouble() - 0.5) * 2 * delta + baseBright;
				var r = Math.Round((c0.R * (1 - t1) + c1.R * t1) * s);
				var g = Math.Round((c0.G * (1 - t1) + c1.G * t1) * s);
				var b = Math.Round((c0.B * (1 - t1) + c1.B * t1) * s);
				r = Math.Min(r, 255);
				g = Math.Min(g, 255);
				b = Math.Min(b, 255);
				var c = Color.FromRgb((byte)r, (byte)g, (byte)b);
				return c;
			}
			return new Color();
		}

		public double CalcPosition(System.Windows.Point p, Geometry2d.Line2d line)
		{
			var pc = new Complex(p.X, p.Y);
			var vec = line.Vector;
			var hp = line.Projection(pc);
			var h = hp - pc;
			var t = -(h.Real + h.Imaginary) / (vec.Real + vec.Imaginary);
			t += (hp - line.begin).Real / vec.Real;
			return t;
		}

		public List<DelaunatorSharp.Point> GeneratePoints(int minX, int maxX, int minY, int maxY, int divX, int divY, int repeatX, int repeatY, int delta)
		{
			var width = maxX - minX;
			var height = maxY - minY;
			var ret = new List<DelaunatorSharp.Point>();
			for (int i = 0; i < divX; i++)
			{
				var x = minX + width / 2.0 / divX * (i * 2 + 1);
				for (int j = 0; j < divY; j++)
				{
					var y = minY + height / 2.0 / divY * (j * 2 + 1);
					var p = new DelaunatorSharp.Point(x + (rand.NextDouble() - 0.5) * 2 * delta, y + (rand.NextDouble() - 0.5) * 2 * delta);
					for (int s = 0; s < repeatX; s++)
					{
						for (int t = 0; t < repeatY; t++)
						{
							ret.Add(new DelaunatorSharp.Point(p.X + width * s, p.Y + height * t));
						}
					}
				}
			}
			return ret;
		}
	}

	// https://qiita.com/tricogimmick/items/894914f6bbe224a45d49
	// Canvas クラスの拡張メソッドとして実装する
	public static class CanvasExtensions
	{
		// Canvas を画像ファイルとして保存する。
		public static void toImage(this Canvas canvas, string path, BitmapEncoder encoder = null)
		{
			// レイアウトを再計算させる
			var size = new Size(canvas.Width, canvas.Height);
			canvas.Measure(size);
			canvas.Arrange(new Rect(size));

			// VisualObjectをBitmapに変換する
			var renderBitmap = new RenderTargetBitmap((int)size.Width,       // 画像の幅
													  (int)size.Height,      // 画像の高さ
													  96.0d,                 // 横96.0DPI
													  96.0d,                 // 縦96.0DPI
													  PixelFormats.Pbgra32); // 32bit(RGBA各8bit)
			renderBitmap.Render(canvas);

			// 出力用の FileStream を作成する
			using (var os = new FileStream(path, FileMode.Create))
			{
				// 変換したBitmapをエンコードしてFileStreamに保存する。
				// BitmapEncoder が指定されなかった場合は、PNG形式とする。
				encoder = encoder ?? new PngBitmapEncoder();
				encoder.Frames.Add(BitmapFrame.Create(renderBitmap));
				encoder.Save(os);
			}
		}
	}

	public static class Geometry2d
	{
		public static double eps = 1e-10;

		/// <summary>
		/// Returns whether <paramref name="x"/> is plus, minus, or zero.
		/// </summary>
		/// <returns>1 if plus, 0 if zero, or -1 if minus.</returns>
		public static int Sign(double x) => (x < eps ? -1 : (x > eps ? 1 : 0));
		public static int Compare(this double x, double y) => Sign(x - y);

		public static int ComplexCompare(this Complex x, Complex y)
		{
			if (Compare(x.Real, y.Real) == 0) return Compare(x.Real, y.Real);
			else return Compare(x.Imaginary, y.Imaginary);
		}

		public static double Distance(this Complex x, Complex y) => (x - y).Magnitude;

		public static Complex Normalize(this Complex x) => x / x.Magnitude;

		public static Complex Rotate(this Complex z, double theta) => z * Complex.FromPolarCoordinates(1, theta);
		public static Complex RotateDegree(this Complex z, double deg) => Rotate(z, deg / 180 * Math.PI);

		public static double Dot(this Complex x, Complex y) => x.Real * y.Real + x.Imaginary * y.Imaginary;
		public static double Cross(this Complex x, Complex y) => x.Real * y.Imaginary - y.Real * x.Imaginary;

		/// <summary>
		/// Returns the relative position of points a, b, c.
		/// </summary>
		/// <returns>
		/// <para>+1 if BC turns left seen from AB.</para>
		/// <para>-1 if BC turns right seen from AB.</para>
		/// <para>+2 if they are on a line in the order ABC or CBA.</para>
		/// <para>0 if they are on a line in the order ACB or BCA.</para>
		/// <para>-2 if they are on a line in the order BAC or CAB.</para></returns>
		public static int ISP(Complex a, Complex b, Complex c)
		{
			var cross = Sign((b - a).Cross(c - b));
			if (cross > 0) return 1;
			else if (cross < 0) return -1;
			else
			{
				var dot = Sign((b - a).Dot(c - a));
				return dot * 2;
			}
		}

		/// <summary>
		/// Returns whether the angle ABC is acute, right, or obtuse.
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		/// <param name="c"></param>
		/// <returns><para>+1 if ABC is an acute angle.</para>
		/// <para>0 if ABC is a right angle.</para>
		/// <para>-1 if ABC is a obtuse angle.</para></returns>
		public static int AngleType(Complex a, Complex b, Complex c) => Sign((a - b).Dot(c - b));

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static double Pow2(double x) => x * x;

		public struct Line2d
		{
			public Complex begin, end;

			/// <summary>
			/// Contructs Line2d that passes the both points.
			/// </summary>
			public Line2d(Complex begin, Complex end)
			{
				this.begin = begin;
				this.end = end;
			}

			/// <summary>
			/// Constructs Line2d using the formula ax+by+c=0.
			/// </summary>
			public Line2d(double a, double b, double c)
			{
				if (Sign(a) == 0 && Sign(b) == 0) throw new ArgumentException();
				if (Sign(a) == 0)
				{
					begin = new Complex(0.0, -c / b);
					end = new Complex(1.0, -c / b);
				}
				else if (Sign(b) == 0)
				{
					begin = new Complex(-c / a, 0.0);
					end = new Complex(-c / a, 1.0);
				}
				else
				{
					begin = new Complex(-c / a, 0.0);
					end = new Complex(0.0, -c / b);
				}
			}

			public Complex Vector { get => end - begin; }
			public Complex CounterVector { get => begin - end; }

			public Complex Intersection(Line2d other)
			{
				if (Sign(this.Vector.Cross(other.Vector)) == 0) return Complex.NaN;
				return this.begin + this.Vector
					* ((other.end - this.begin).Cross(other.Vector) / this.Vector.Cross(other.Vector));
			}

			public double Distance(Complex other) => Math.Abs(this.Vector.Cross(other - this.begin) / this.Vector.Magnitude);

			public Complex Projection(Complex other)
				=> this.begin + this.Vector.Normalize() * this.Vector.Dot(other - this.begin) / this.Vector.Magnitude;

			public Complex Reflection(Complex other) => other + 2 * (this.Projection(other) - other);
		}

		public class Circle2d
		{
			Complex center;
			double radius;

			public Circle2d(Complex center, double radius)
			{
				this.center = center;
				this.radius = radius;
			}

			public Complex[] Intersection(Line2d other)
			{
				var comp = other.Distance(center).Compare(radius);
				if (comp < 0) return null;
				else if (comp == 0)
				{
					var ret = new Complex[1];
					ret[0] = other.Projection(center);
					return ret;
				}
				else
				{
					var ret = new Complex[2];
					var h = other.Projection(center);
					var ohLen = (h - center).Magnitude;
					var vec = (h - center).RotateDegree(90);
					var len = Math.Sqrt(Math.Max(0, radius * radius - ohLen * ohLen));
					ret[0] = h + vec.Normalize() * len;
					ret[1] = h - vec.Normalize() * len;
					return ret;
				}
			}

			public Complex[] Intersection(Circle2d other)
			{
				var dist = (this.center).Distance(other.center);
				var comp = dist.Compare(this.radius + other.radius);
				if (comp > 0) return null;
				else
				{
					var x = (Pow2(dist) + Pow2(this.radius) - Pow2(other.radius)) / (2 * dist);
					var h = (other.center - this.center).Normalize() * x;
					if (comp == 1)
					{
						var ret = new Complex[1];
						ret[0] = h;
						return ret;
					}
					var vec = (other.center - this.center).RotateDegree(90);
					var y = Math.Sqrt(Math.Max(0, Pow2(this.radius) - Pow2(x)));
					var ret2 = new Complex[2];
					ret2[0] = h + vec.Normalize() * y;
					ret2[1] = h - vec.Normalize() * y;
					return ret2;
				}
			}
		}
	}
}
