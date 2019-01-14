class Planet
{
public:
  vec mass;
  vec initial_position;
  vec intitial_velocity;
  vec position;
  vec velocity;

  Planet(m, x0, y0, z0, vx0, vy0, vz0);

  void change_initial_position(x_new, y_new, z_new);
  void change_initial_velocity(vx_new, vy_new, vz_new);
  void change_mass(m_new);

};
