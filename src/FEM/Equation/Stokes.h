//*****************************************************************************
//Title		:src/FEM/Equation/Stokes.h
//Author	:Tanabe Yuta
//Date		:2020/02/05
//Copyright	:(C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>


#include "../../LinearAlgebra/Models/Matrix.h"
#include "../../LinearAlgebra/Models/Vector.h"


namespace PANSFEM2 {
    //******************************Get element stiffness matrix******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void Stokes(Matrix<T>& _Ke, std::vector<Vector<T> >& _x, std::vector<int>& _elementu, std::vector<int>& _elementp, T _mu) {
		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p
        
        //----------Initialize element matrix----------
		_Ke = Matrix<T>(2*m + n, 2*m + n);

		//----------Generate cordinate matrix X of velosity u----------
		Matrix<T> Xu = Matrix<T>(0, 2);
		for(int i = 0; i < m; i++){
			Xu = Xu.Vstack(_x[_elementu[i]].Transpose());
		}

        //----------Generate cordinate matrix X of pressure p----------
		Matrix<T> Xp = Matrix<T>(0, 2);
		for(int i = 0; i < n; i++){
			Xp = Xp.Vstack(_x[_elementp[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
            //----------Get shape function for pressure p----------
			Vector<T> N = SFP<T>::N(IC<T>::Points[g]);
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;

			//----------Get shape function for velocity u----------
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			Matrix<T> dMdr = SFU<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dMdX = dXdr.Inverse()*dMdr;

            //----------Get K matrix----------
            Matrix<T> K = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    K(3*i, 3*j) = 2.0*_mu*dMdX(0, i)*dMdX(0, j) + _mu*dMdX(1, i)*dMdX(1, j);    K(3*i, 3*j + 1) = _mu*dMdX(0, i)*dMdX(1, j);                                        K(3*i, 3*j + 2) = -N(j)*dMdX(0, i);
                    K(3*i + 1, 3*j) = _mu*dMdX(1, i)*dMdX(0, j);                                K(3*i + 1, 3*j + 1) = _mu*dMdX(0, i)*dMdX(0, j) + 2.0*_mu*dMdX(1, i)*dMdX(1, j);    K(3*i + 1, 3*j + 2) = -N(j)*dMdX(1, i);
                    K(3*i + 2, 3*j) = -N(i)*dMdX(0, j);                                         K(3*i + 2, 3*j + 1) = -N(i)*dMdX(1, j);                                             K(3*i + 2, 3*j + 2) = T();
                }
                for(int j = n; j < m; j++){
                    K(3*i, n + 2*j) = 2.0*_mu*dMdX(0, i)*dMdX(0, j) + _mu*dMdX(1, i)*dMdX(1, j);    K(3*i, n + 2*j + 1) = _mu*dMdX(0, i)*dMdX(1, j);
                    K(3*i + 1, n + 2*j) = _mu*dMdX(1, i)*dMdX(0, j);                                K(3*i + 1, n + 2*j + 1) = _mu*dMdX(0, i)*dMdX(0, j) + 2.0*_mu*dMdX(1, i)*dMdX(1, j);
                    K(3*i + 2, n + 2*j) = -N(i)*dMdX(0, j);                                         K(3*i + 2, n + 2*j + 1) = -N(i)*dMdX(1, j);
                }
            }
            for(int i = n; i < m; i++){
                for(int j = 0; j < n; j++){
                    K(n + 2*i, 3*j) = 2.0*_mu*dMdX(0, i)*dMdX(0, j) + _mu*dMdX(1, i)*dMdX(1, j);    K(n + 2*i, 3*j + 1) = _mu*dMdX(0, i)*dMdX(1, j);                                        K(n + 2*i, 3*j + 2) = -N(j)*dMdX(0, i);
                    K(n + 2*i + 1, 3*j) = _mu*dMdX(1, i)*dMdX(0, j);                                K(n + 2*i + 1, 3*j + 1) = _mu*dMdX(0, i)*dMdX(0, j) + 2.0*_mu*dMdX(1, i)*dMdX(1, j);    K(n + 2*i + 1, 3*j + 2) = -N(j)*dMdX(1, i);
                }
                for(int j = n; j < m; j++){
                    K(n + 2*i, n + 2*j) = 2.0*_mu*dMdX(0, i)*dMdX(0, j) + _mu*dMdX(1, i)*dMdX(1, j);    K(n + 2*i, n + 2*j + 1) = _mu*dMdX(0, i)*dMdX(1, j);
                    K(n + 2*i + 1, n + 2*j) = _mu*dMdX(1, i)*dMdX(0, j);                                K(n + 2*i + 1, n + 2*j + 1) = _mu*dMdX(0, i)*dMdX(0, j) + 2.0*_mu*dMdX(1, i)*dMdX(1, j);
                }
            }
			
			//----------Update element matrix----------
			_Ke += K*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}


    //******************************Get element traction vector******************************
    template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
    void StokesTraction(Vector<T>& _Fe, std::vector<Vector<T> >& _x, std::vector<int>& _elementu, std::vector<int>& _elementp, T _sx, T _sy){
        int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p
        
        //----------Initialize element matrix----------
		_Fe = Vector<T>(2*m + n);

		//----------Generate cordinate matrix X of velosity u----------
		Matrix<T> Xu = Matrix<T>(0, 2);
		for(int i = 0; i < m; i++){
			Xu = Xu.Vstack(_x[_elementu[i]].Transpose());
		}

        //----------Generate cordinate matrix X of pressure p----------
		Matrix<T> Xp = Matrix<T>(0, 2);
		for(int i = 0; i < n; i++){
			Xp = Xp.Vstack(_x[_elementp[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
            //----------Get shape function for pressure p----------
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
            T dl = sqrt((dXdr*dXdr.Transpose())(0, 0));

			//----------Get shape function for velocity u----------
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);
			
            //----------Geenerate B matrix----------
            Matrix<T> B = Matrix<T>(2, 2*m + n);
            for(int i = 0; i < n; i++){
                B(0, 3*i) = M(i);   B(0, 3*i + 1) = T();    B(0, 3*i + 2) = T();
                B(1, 3*i) = T();    B(1, 3*i + 1) = M(i);   B(1, 3*i + 2) = T();
            }
            for(int i = n; i < m; i++){
                B(0, n + 2*i) = M(i);   B(0, n + 2*i + 1) = T();
                B(1, n + 2*i) = T();    B(1, n + 2*i + 1) = M(i);
            }
            
            //----------Get surface force vector----------
			Vector<T> s = Vector<T>({ _sx, _sy });
			
			//----------Update element matrix----------
			_Fe += B.Transpose()*s*dl*IC<T>::Weights[g][0];
		}
    }


    //******************************Get element mass matrix******************************
	template<class T, template<class>class SFU, template<class>class SFP, template<class>class IC>
	void StokesMass(Matrix<T>& _Ce, std::vector<Vector<T> >& _x, std::vector<int>& _elementu, std::vector<int>& _elementp, T _rho) {
		int m = _elementu.size();   //  Number of shapefunction for velosity u
        int n = _elementp.size();   //  Number of shapefunction for pressure p
        
        //----------Initialize element matrix----------
		_Ce = Matrix<T>(2*m + n, 2*m + n);

		//----------Generate cordinate matrix X of velosity u----------
		Matrix<T> Xu = Matrix<T>(0, 2);
		for(int i = 0; i < m; i++){
			Xu = Xu.Vstack(_x[_elementu[i]].Transpose());
		}

        //----------Generate cordinate matrix X of pressure p----------
		Matrix<T> Xp = Matrix<T>(0, 2);
		for(int i = 0; i < n; i++){
			Xp = Xp.Vstack(_x[_elementp[i]].Transpose());
		}

		//----------Loop of Gauss Integration----------
		for (int g = 0; g < IC<T>::N; g++) {
            //----------Get shape function for pressure p----------
			Matrix<T> dNdr = SFP<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*Xp;
			T J = dXdr.Determinant();
			
			//----------Get shape function for velocity u----------
			Vector<T> M = SFU<T>::N(IC<T>::Points[g]);

            //----------Get K matrix----------
            Matrix<T> C = Matrix<T>(2*m + n, 2*m + n);
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    C(3*i, 3*j) = _rho*M(i)*M(j);   C(3*i, 3*j + 1) = T();                  C(3*i, 3*j + 2) = T();
                    C(3*i + 1, 3*j) = T();          C(3*i + 1, 3*j + 1) = _rho*M(i)*M(j);   C(3*i + 1, 3*j + 2) = T();
                    C(3*i + 2, 3*j) = T();          C(3*i + 2, 3*j + 1) = T();              C(3*i + 2, 3*j + 2) = T();
                }
                for(int j = n; j < m; j++){
                    C(3*i, n + 2*j) = _rho*M(i)*M(j);   C(3*i, n + 2*j + 1) = T();
                    C(3*i + 1, n + 2*j) = T();          C(3*i + 1, n + 2*j + 1) = _rho*M(i)*M(j);
                }
            }
            for(int i = n; i < m; i++){
                for(int j = 0; j < n; j++){
                    C(n + 2*i, 3*j) = _rho*M(i)*M(j);   C(n + 2*i, 3*j + 1) = T();
                    C(n + 2*i + 1, 3*j) = T();          C(n + 2*i + 1, 3*j + 1) = _rho*M(i)*M(j);
                }
                for(int j = n; j < m; j++){
                    C(n + 2*i, n + 2*j) = _rho*M(i)*M(j);   C(n + 2*i, n + 2*j + 1) = T();
                    C(n + 2*i + 1, n + 2*j) = T();          C(n + 2*i + 1, n + 2*j + 1) = _rho*M(i)*M(j);
                }
            }
			
			//----------Update element matrix----------
			_Ce += C*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
		}
	}
}