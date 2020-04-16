#include "fem.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>

namespace FEM
{

DeformableMesh2D::DeformableMesh2D(std::vector<float> nodes_x, std::vector<float> nodes_y, std::vector<Element2D> elements, float poisson_ratio, float young_modulus) {
	D_ <<
		1.0f,           poisson_ratio, 0.0f,
		poisson_ratio, 1.0f,           0.0f,
		0.0f,           0.0f,           (1.0f - poisson_ratio) / 2.0f;

	D_ *= young_modulus / (1.0f - poisson_ratio * poisson_ratio);

	nodes_x_ = nodes_x;
	nodes_y_ = nodes_y;

	elements_ = elements;

	nodes_count_ = nodes_x.size();

	forces_.resize(2 * nodes_count_);
	forces_.setZero();

	Q_.resize(2 * nodes_count_, 2 * nodes_count_);
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

	Q_.setFromTriplets(triplets.begin(), triplets.end());

	// TODO: Apply constraints.
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

	for (int k = 0; k < Q_.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<float>::InnerIterator it(Q_, k); it; ++it)
		{
			for (std::vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit)
			{
				SetConstraints(it, *idit);
			}
		}
	}

	std::cout << "Global stiffness matrix:\n";
	std::cout << Q_ << std::endl;

	// Invert Q to get compliance matrix
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver(Q_);
	Eigen::SparseMatrix<float> I(2 * nodes_count_, 2 * nodes_count_);
	I.setIdentity();
	Q_ = solver.solve(I);
}

void DeformableMesh2D::setForce(int node, float x, float y) {
	forces_[2 * node + 0] = x;
	forces_[2 * node + 1] = y;
};

Eigen::VectorXf DeformableMesh2D::calculateDisplacements() {
	Eigen::VectorXf displacements = Q_ * forces_;
	return displacements;
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