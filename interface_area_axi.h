double interface_area_axi (scalar c)
 {
  double area = 0.;
  foreach (reduction(+:area))
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = interface_normal (point, c), p;
      double alpha = plane_alpha (c[], n);
      area += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p) *2*M_PI*y;
    }
  return area;
}