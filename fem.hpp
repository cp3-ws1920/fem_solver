#include <array>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

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
 private:
  int nodes_count_;

  std::vector<Element2D> elements_;
  std::vector<float> nodes_x_;
  std::vector<float> nodes_y_;

  Eigen::VectorXf forces_;
  
  std::vector<Constraint> constraints_;

  Eigen::Matrix3f D_; // Elasticity matrix.
  Eigen::SparseMatrix<float> Q_; // Compliance matrix (inverse of stiffness matrix.)

  void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index);
 public:
  DeformableMesh2D(std::vector<float> nodes_x, std::vector<float> nodes_y, std::vector<Element2D> elements, float poisson_ratio, float young_modulus);
  // Pre-processing.
  void setConstraint(Constraint constraint);
  void calculateMatrix();
  // Getting displacements from forces.
  void setForce(int node, float x, float y);
  Eigen::VectorXf calculateDisplacements();
};

}