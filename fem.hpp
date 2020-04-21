#ifndef FEM_HPP_
#define FEM_HPP_

#include <array>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace FEM
{

class Element2D {
 private:
  Eigen::Matrix<float, 3, 6> B_;
 public:
  std::array<int, 3> node_ids_; // TODO: This should be private
  void CalculateStiffnessMatrix(const Eigen::Matrix3f& D, const std::vector<float>& nodes_x,
                                const std::vector<float>& nodes_y, std::vector<Eigen::Triplet<float> >& triplets);
};

struct Constraint
{
  enum Type
  {
    NONE = 0,
    UX = 1 << 0,
    UY = 1 << 1,
    UXY = UX | UY
  };
  Type type;
  int node;
};

class DeformableMesh2D {
 public:
  enum Method {
		METHOD_EXPLICIT_EULER,
		METHOD_MODIFIED_EXPLICIT_EULER,
		METHOD_IMPROVED_EULER,
    METHOD_RUNGE_KUTTA_3,
    METHOD_RUNGE_KUTTA_4,
    METHOD_IMPLICIT_EULER
	};
 private:
  int nodes_count_;

  float poisson_ratio_;
  float young_modulus_;

  std::vector<Element2D> elements_;
  std::vector<float> nodes_x_;
  std::vector<float> nodes_y_;

  Method method_;
  float fixed_delta_;
  bool fixed_delta_enabled_;

  Eigen::VectorXf forces_;
  
  std::vector<Constraint> constraints_;

  Eigen::Matrix3f D_; // Elasticity matrix.
  Eigen::SparseMatrix<float> Q_; // Compliance matrix (inverse of stiffness matrix.)

  // Used for motion equations
  Eigen::SparseMatrix<float> K_; // Stiffnes matrix
  Eigen::VectorXf D1_; // Displacements
  Eigen::VectorXf D2_; // Velocities
  float R_; // Friction

  // Procedure matrices for implicit euler with fixed delta
  Eigen::SparseMatrix<float> T1_, T2_, T3_, T4_;

  void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index);
  void calculateMatrix();
  void calculateProcedureMatrix();
 public:
  DeformableMesh2D();
  void preprocess();
  // Pre-processing.
  void clearContraints() { constraints_.clear(); }
  void setConstraint(Constraint constraint);
  // Getting displacements from forces.
  void setYoungModulus(float young_modulus) { young_modulus_ = young_modulus; }
  void setPoissonRatio(float poisson_ratio) { poisson_ratio_ = poisson_ratio; }
  void setForce(int node, float x, float y);
  void setFriction(float R) { R_ = R; }
  float getFriction() { return R_; }
  std::vector<float> calculateDisplacements();
  std::vector<float> freeOscillationStep(float delta);
  void resetVelocity();
  void setElements(std::vector<Element2D> elements) { elements_ = elements; }
  std::vector<Element2D> getElements() { return elements_; }
  void setNodesX(std::vector<float> nodes_x) { nodes_x_ = nodes_x; }
  std::vector<float> getNodesX() { return nodes_x_; }
  void setNodesY(std::vector<float> nodes_y) { nodes_y_ = nodes_y; }
  std::vector<float> getNodesY() { return nodes_y_; }
  std::vector<float> getDisplacements() { return std::vector<float>(D1_.data(), D1_.data() + D1_.size()); }
  std::vector<float> getVelocities() { return std::vector<float>(D2_.data(), D2_.data() + D2_.size()); };
  void setMethod(Method method) { method_ = method; }
  Method getMethod() { return method_; }
  void setFixedDeltaEnabled(bool enabled) { fixed_delta_enabled_ = enabled; }
  bool getFixedDeltaEnabled() { return fixed_delta_enabled_; }
  void setFixedDelta(float fixed_delta) { fixed_delta_ = fixed_delta; }
  float getFixedDelta() { return fixed_delta_; }
};

}

#endif // FEM_HPP_