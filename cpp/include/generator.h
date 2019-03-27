//
// Created by carlos on 19/03/19.
//

#ifndef ONLINE_PLANNING_GENERATOR_H
#define ONLINE_PLANNING_GENERATOR_H

#include "bezier.h"
#include "model.h"

struct PhysLimits {
    Eigen::VectorXd pmax;
    Eigen::VectorXd pmin;
    Eigen::VectorXd amax;
    Eigen::VectorXd amin;
};

struct TuningParams {
    float s_free, s_obs, s_repel;
    int spd_f, spd_o, spd_r;
    double lin_coll, quad_coll;
    Eigen::VectorXd energy_weights;
};

struct CollisionParams {
    int order;
    float rmin;
    float c;
};

struct MpcParams {
    float h, Ts;
    int k_hor;
    const TuningParams& tuning;
    const PhysLimits& limits;
    const CollisionParams coll;
};

class Generator {
public:
    struct Params {
        const BezierCurve::Params& bezier_params;
        const DoubleIntegrator3D::Params& model_params;
        const MpcParams& mpc_params;
        const Eigen::MatrixXd& po;
        const Eigen::MatrixXd& pf;
    };

    Generator(const Generator::Params& p);
    ~Generator(){};

private:
    float _h;
    float _Ts;
    float _k_hor;
    int _l;
    int _num_ctrl_pts;
    int _dim;

    Eigen::MatrixXd _po;
    Eigen::MatrixXd _pf;

    BezierCurve _bezier;
    DoubleIntegrator3D _model_pred;
    DoubleIntegrator3D _model_exec;
    Eigen::MatrixXd _H_energy;
    Constraint _ineq;
    Constraint _eq;
    StatePropagator _Phi_pred;
    StatePropagator _Phi_exec;

    // Matrices to minimize goal error
    Eigen::MatrixXd _H_f;
    Eigen::MatrixXd _H_o;
    Eigen::MatrixXd _H_r;

    // Methods
    Constraint build_ineq_constr(const PhysLimits& lim);
    void set_error_penalty_mats(const TuningParams& p, const Eigen::MatrixXd& pf);
};

#endif //ONLINE_PLANNING_GENERATOR_H