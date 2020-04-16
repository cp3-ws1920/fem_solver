#include "fem.hpp"

#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace FEM
{

DeformableMesh2D::DeformableMesh2D() {
	fixed_delta_ = 0.0f; 
	fixed_delta_enabled_ = false;
}

void DeformableMesh2D::preprocess() {
	D_ <<
		1.0f,           poisson_ratio_, 0.0f,
		poisson_ratio_, 1.0f,           0.0f,
		0.0f,           0.0f,           (1.0f - poisson_ratio_) / 2.0f;

	D_ *= young_modulus_ / (1.0f - poisson_ratio_ * poisson_ratio_);

	nodes_count_ = nodes_x_.size();

	forces_.resize(2 * nodes_count_);
	forces_.setZero();

	Q_.resize(2 * nodes_count_, 2 * nodes_count_);
	K_.resize(2 * nodes_count_, 2 * nodes_count_);

	D2_.resize(2 * nodes_count_);
	D1_.resize(2 * nodes_count_);
	D1_.setZero();
	D2_.setZero();

	calculateMatrix();
	calculateProcedureMatrix();
}

void DeformableMesh2D::calculateProcedureMatrix() {
	Eigen::SparseMatrix<float> I(2 * nodes_count_, 2 * nodes_count_);
	I.setIdentity();
	Eigen::SparseMatrix<float> S = I * (1 + fixed_delta_ * R_) - fixed_delta_*fixed_delta_ * K_;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver(S);
	Eigen::SparseMatrix<float> S_inv = solver.solve(I);
	T1_ = (I + fixed_delta_*fixed_delta_ * K_ * S_inv);
	T2_ = fixed_delta_ * S_inv;
	T3_ = fixed_delta_ * K_ * S_inv;
	T4_ = S_inv;
}

void DeformableMesh2D::setConstraint(Constraint constraint) {
	constraints_.push_back(constraint);
}

void DeformableMesh2D::SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index)
{
	if (it.row() == index || it.col() == index)
	{
		it.valueRef() = it.row() == it.col() ? 1.0f : 0.0f;
	}
};

void DeformableMesh2D::calculateMatrix()
{
	std::vector<Eigen::Triplet<float> > triplets;
	for (std::vector<Element2D>::iterator it = elements_.begin(); it != elements_.end(); ++it) {
		it->CalculateStiffnessMatrix(D_, nodes_x_, nodes_y_, triplets);
	}

	K_.setFromTriplets(triplets.begin(), triplets.end());

	std::vector<int> indicesToConstraint;

	for (std::vector<Constraint>::const_iterator it = constraints_.begin(); it != constraints_.end(); ++it)
	{
		if (it->type & Constraint::UX)
		{
			indicesToConstraint.push_back(2 * it->node + 0);
		}
		if (it->type & Constraint::UY)
		{
			indicesToConstraint.push_back(2 * it->node + 1);
		}
	}

	for (int k = 0; k < K_.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<float>::InnerIterator it(K_, k); it; ++it)
		{
			for (std::vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit)
			{
				SetConstraints(it, *idit);
			}
		}
	}

	// Invert Q to get compliance matrix
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver(K_);
	Eigen::SparseMatrix<float> I(2 * nodes_count_, 2 * nodes_count_);
	I.setIdentity();
	Q_ = solver.solve(I);
}

void DeformableMesh2D::setForce(int node, float x, float y) {
	forces_[2 * node + 0] = x;
	forces_[2 * node + 1] = y;
};

Eigen::VectorXf DeformableMesh2D::calculateDisplacements() {
	D1_ = Q_ * forces_;
	return D1_;
}

void DeformableMesh2D::resetVelocity() {
	D2_.setZero();
}

Eigen::VectorXf DeformableMesh2D::freeOcillationStep(float delta = 0.0f) {
	if (fixed_delta_enabled_) delta = fixed_delta_;
	switch (method_) {
		case (DeformableMesh2D::Method::METHOD_EXPLICIT_EULER) : {
			Eigen::VectorXf k1_1 = D2_;
			Eigen::VectorXf k1_2 = K_ * D1_ - R_ * D2_;
			D1_ = D1_ + delta * k1_1;
			D2_ = D2_ + delta * k1_2;
			break;
		}
		case (DeformableMesh2D::Method::METHOD_MODIFIED_EXPLICIT_EULER) : {
			D2_ = D2_ + delta * K_ * D1_ - delta * R_ * D2_;
			D1_ = D1_ + delta * D2_;
			break;
		}
		case (DeformableMesh2D::Method::METHOD_IMPROVED_EULER) : {
			Eigen::VectorXf k1_1 = D2_;
			Eigen::VectorXf k1_2 = K_ * D1_ - R_ * D2_;
			Eigen::VectorXf k2_1 = D2_ + delta / 2.0f * k1_2;
			Eigen::VectorXf k2_2 = K_ * (D1_ + delta / 2.0f * k1_1) - R_ * (D2_ + delta / 2.0f * k1_2);
			D1_ = D1_ + delta * k1_1;
			D2_ = D2_ + delta * k1_2;
			break;
		}
		case (DeformableMesh2D::Method::METHOD_RUNGE_KUTTA_3) : {
			Eigen::VectorXf k1_1 = D2_;
			Eigen::VectorXf k1_2 = K_ * D1_ - R_ * D2_;
			Eigen::VectorXf k2_1 = D2_ + delta / 2.0f * k1_2;
			Eigen::VectorXf k2_2 = K_ * (D1_ + delta / 2.0f * k1_1) - R_ * (D2_ + delta / 2.0f * k1_2);
			Eigen::VectorXf k3_1 = D2_ - delta * k1_2 + 2.0f * delta * k2_2;
			Eigen::VectorXf k3_2 = K_ * (D1_ - delta * k1_1 + 2.0f * delta * k2_1) - R_ * (D2_ - delta * k1_2 +  2.0f * delta * k2_2);
			D1_ = D1_ + delta * (k1_1 + 4.0f * k2_1 + k3_1) / 6.0f;
			D2_ = D2_ + delta * (k1_2 + 4.0f * k2_2 + k3_2) / 6.0f;
			break;
		}
		case (DeformableMesh2D::Method::METHOD_RUNGE_KUTTA_4) : {
			Eigen::VectorXf k1_1 = D2_;
			Eigen::VectorXf k1_2 = K_ * D1_ - R_ * D2_;
			Eigen::VectorXf k2_1 = D2_ + delta / 2.0f * k1_2;
			Eigen::VectorXf k2_2 = K_ * (D1_ + delta / 2.0f * k1_1) - R_ * (D2_ + delta / 2.0f * k1_2);
			Eigen::VectorXf k3_1 = D2_ + delta / 2.0f * k2_2;
			Eigen::VectorXf k3_2 = K_ * (D1_ + delta / 2.0f * k2_1) - R_ * (D2_ + delta / 2.0f * k2_2);
			Eigen::VectorXf k4_1 = D2_ + delta * k2_2;
			Eigen::VectorXf k4_2 = K_ * (D1_ + delta * k3_1) - R_ * (D2_ + delta * k3_2);
			D1_ = D1_ + delta * (k1_1 + 2.0f * k2_1 + 2.0f * k3_1 + k4_1) / 6.0f;
			D2_ = D2_ + delta * (k1_2 + 2.0f * k2_2 + 2.0f * k3_2 + k4_2) / 6.0f;
			break;
		}
		case (DeformableMesh2D::Method::METHOD_IMPLICIT_EULER) : {
			Eigen::VectorXf D1 = D1_;
			Eigen::VectorXf D2 = D2_;
			if (!fixed_delta_enabled_)
			{
				Eigen::SparseMatrix<float> I(2 * nodes_count_, 2 * nodes_count_);
				I.setIdentity();
				Eigen::SparseMatrix<float> S = I * (1 + delta * R_) - delta * delta * K_;
				Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver(S);
				Eigen::SparseMatrix<float> S_inv = solver.solve(I);
				D1_ = (I + delta*delta * K_ * S_inv) * D1 + delta * S_inv * D2;
				D2_ =  delta * K_ * S_inv * D1 + S_inv * D2;
			} else {
				D1_ = T1_ * D1 + T2_ * D2;
				D2_ = T3_ * D1 + T4_ * D2;
			}
			break;
		}
	}

	return D1_;
}

void Element2D::CalculateStiffnessMatrix(const Eigen::Matrix3f& D, const std::vector<float>& nodes_x,
                                       const std::vector<float>& nodes_y, std::vector<Eigen::Triplet<float> >& triplets)
{
	Eigen::Vector3f x, y;
	x << nodes_x[node_ids_[0]], nodes_x[node_ids_[1]], nodes_x[node_ids_[2]];
	y << nodes_y[node_ids_[0]], nodes_y[node_ids_[1]], nodes_y[node_ids_[2]];

	
	Eigen::Matrix3f C;
	C << Eigen::Vector3f(1.0f, 1.0f, 1.0f), x, y;
	
	Eigen::Matrix3f IC = C.inverse();

	for (int i = 0; i < 3; i++)
	{
		B_(0, 2 * i + 0) = IC(1, i);
		B_(0, 2 * i + 1) = 0.0f;
		B_(1, 2 * i + 0) = 0.0f;
		B_(1, 2 * i + 1) = IC(2, i);
		B_(2, 2 * i + 0) = IC(2, i);
		B_(2, 2 * i + 1) = IC(1, i);
	}
	Eigen::Matrix<float, 6, 6> K = B_.transpose() * D * B_ * C.determinant() / 2.0f;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Eigen::Triplet<float> trplt11(2 * node_ids_[i] + 0, 2 * node_ids_[j] + 0, K(2 * i + 0, 2 * j + 0));
			Eigen::Triplet<float> trplt12(2 * node_ids_[i] + 0, 2 * node_ids_[j] + 1, K(2 * i + 0, 2 * j + 1));
			Eigen::Triplet<float> trplt21(2 * node_ids_[i] + 1, 2 * node_ids_[j] + 0, K(2 * i + 1, 2 * j + 0));
			Eigen::Triplet<float> trplt22(2 * node_ids_[i] + 1, 2 * node_ids_[j] + 1, K(2 * i + 1, 2 * j + 1));

			triplets.push_back(trplt11);
			triplets.push_back(trplt12);
			triplets.push_back(trplt21);
			triplets.push_back(trplt22);
		}
	}
}

}