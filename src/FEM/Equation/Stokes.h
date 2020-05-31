//*****************************************************************************
//  Title		:   src/FEM/Equation/Stokes.h
//  Author	    :   Tanabe Yuta
//  Date		:   2020/02/17
//  Copyright	:   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //******************************Get element stiffness matrix******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void StokesStiffness(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _mu) {
        assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xu = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			Xu(i, 0) = _x[_elementu[i]](0); Xu(i, 1) = _x[_elementu[i]](1);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    K(i, j) = 2.0*_mu*dMdX(0, i)*dMdX(0, j) + _mu*dMdX(1, i)*dMdX(1, j);    K(i, j + m) = _mu*dMdX(0, i)*dMdX(1, j);
                    K(i + m, j) = _mu*dMdX(1, i)*dMdX(0, j);                                K(i + m, j + m) = _mu*dMdX(0, i)*dMdX(0, j) + 2.0*_mu*dMdX(1, i)*dMdX(1, j);  
                }
                for(int j = 0; j < n; j++) {
                    K(i, j + 2*m) = -N(j)*dMdX(0, i);
                    K(i + m, j + 2*m) = -N(j)*dMdX(1, i);
                    K(j + 2*m, i) = -N(j)*dMdX(0, i);
                    K(j + 2*m, i + m) = -N(j)*dMdX(1, i);
                }
            }

			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element stiffness matrix of Stokes equation with Brinkman model******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void BrinkmanStokesStiffness(Matrix<T>& _Ke, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _mu, Matrix<T> _k) {
        assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ke = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xu = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			Xu(i, 0) = _x[_elementu[i]](0); Xu(i, 1) = _x[_elementu[i]](1);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < m; j++) {
                    K(i, j) = 2.0*_mu*dMdX(0, i)*dMdX(0, j) + _mu*dMdX(1, i)*dMdX(1, j) + M(i)*_k(0, 0)*M(j);	K(i, j + m) = _mu*dMdX(0, i)*dMdX(1, j) + M(i)*_k(0, 1)*M(j);
                    K(i + m, j) = _mu*dMdX(1, i)*dMdX(0, j) + M(i)*_k(1, 0)*M(j);                             	K(i + m, j + m) = _mu*dMdX(0, i)*dMdX(0, j) + 2.0*_mu*dMdX(1, i)*dMdX(1, j) + M(i)*_k(1, 1)*M(j);  
                }
                for(int j = 0; j < n; j++) {
                    K(i, j + 2*m) = -N(j)*dMdX(0, i);
                    K(i + m, j + 2*m) = -N(j)*dMdX(1, i);
                    K(j + 2*m, i) = -N(j)*dMdX(0, i);
                    K(j + 2*m, i + m) = -N(j)*dMdX(1, i);
                }
            }

			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


    //******************************Get element mass matrix******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void StokesMass(Matrix<T>& _Ce, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, T _rho) {
		assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Ce = Matrix<T>(2*m + n, 2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

		 Matrix<T> Xu = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			Xu(i, 0) = _x[_elementu[i]](0); Xu(i, 1) = _x[_elementu[i]](1);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);

            Matrix<T> C = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    C(i, j) = _rho*M(i)*M(j);   C(i, j + m) = T();
                    C(i + m, j) = T();          C(i + m, j + m) = _rho*M(i)*M(j);
                }
            }
			
			_Ce += C*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


	//******************************Get element traction vector******************************
    template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC, class F>
    void StokesSurfaceForce(Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, F _f){
        assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Fe = Vector<T>(2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xu = Matrix<T>(m, 2);
		for(int i = 0; i < m; i++){
			Xu(i, 0) = _x[_elementu[i]](0); Xu(i, 1) = _x[_elementu[i]](1);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
            Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
            Vector<T> x = Xp.Transpose()*N;
			Matrix<T> dXdr = dNdr*Xp;
            T dl = sqrt((dXdr*dXdr.Transpose())(0, 0));

			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
		
            Matrix<T> B = Matrix<T>(2, 2*m + n);
            for(int i = 0; i < m; i++){
                B(0, i) = M(i);   B(0, i + m) = T();
                B(1, i) = T();    B(1, i + m) = M(i);
            }
            
			_Fe += B.Transpose()*_f(x)*dl*IC<T>::Weights[g][0];
		}
    }


	//******************************Get element body force vector******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC, class F>
	void StokesBodyForce(Vector<T>& _Fe, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementu, const std::vector<int>& _elementu, std::vector<std::vector<std::pair<int, int> > >& _nodetoelementp, const std::vector<int>& _elementp, const std::vector<int>& _doulist, std::vector<Vector<T> >& _x, F _f) {
        assert(_doulist.size() == 3);

		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p

		_Fe = Vector<T>(2*m + n);
		_nodetoelementu = std::vector<std::vector<std::pair<int, int> > >(m, std::vector<std::pair<int, int> >(2));
		for(int i = 0; i < m; i++) {
			_nodetoelementu[i][0] = std::make_pair(_doulist[0], i);
			_nodetoelementu[i][1] = std::make_pair(_doulist[1], m + i);
		}
        _nodetoelementp = std::vector<std::vector<std::pair<int, int> > >(n, std::vector<std::pair<int, int> >(1));
		for(int i = 0; i < n; i++) {
			_nodetoelementp[i][0] = std::make_pair(_doulist[2], 2*m + i);
		}

        Matrix<T> Xp = Matrix<T>(n, 2);
		for(int i = 0; i < n; i++){
			Xp(i, 0) = _x[_elementp[i]](0); Xp(i, 1) = _x[_elementp[i]](1);
		}

		for (int g = 0; g < IC<T>::N; g++) {
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Vector<T> x = Xp.Transpose()*N;
			Vector<T> b = _f(x);

            Vector<T> Fe = Vector<T>(2*m + n);
            for(int i = 0; i < m; i++) {
                Fe(2*i) = b(0)*M(i);	Fe(2*i + 1) = b(1)*M(i);
            }

			_Fe += Fe*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}