/*!
 * Compute the winding number of a boundary with respect to a given coordinate.
 */

#pragma once

#include <cmath>

// FIXME: Want to test T's attributes with type traits to allow T to be cv::Point or
// morph::Vector or morph::vVector or std::array etc etc.
//
#include <type_traits>

namespace morph {

    /*!
     * A winding number class
     *
     * \tparam T the (2D) coordinate type (this might be cv::Point)
     *
     * \tparam Container Something like an std::vector, std::list or std::array,
     * containing a path of points.
     */
    template < typename T,
               template <typename, typename> typename Container,
               typename TT=T,
               typename Alloc=std::allocator<TT> >
    class Winder
    {
    public:
        //! Construct with the boundary reference.
        Winder (const Container<T, Alloc>& _boundary)
            : boundary(_boundary) {}


        //! Compute the winding number of the coordinate px with respect to the boundary.
        int wind (const T& px) {
            this->reset();
            for (auto bp : this->boundary) {
                this->wind (px, bp);
            }
            // Do first pixel again to complete the winding:
            T firstpoint = this->boundary.front();
            this->wind (px, firstpoint);
            double winding_no_d = std::round (this->angle_sum/morph::TWO_PI_D);
            int winding_no = static_cast<int>(winding_no_d);
            return winding_no;
        }

    private:
        //template <typename _T=T, std::enable_if_t<std::is_object<std::decay_t<typename _T::x>>::value, int> = 1 >
        void wind (const T& px, const T& bp) {

            // Get angle from px to bp.

            // This code assumes that T can support operator-
            T pt = bp - px;

            double raw_angle = 0.0;
            // Could try compile-time if here? https://www.bfilipek.com/2016/02/notes-on-c-sfinae.html
            if constexpr (cv::traits::Type<T>::value > -1) { // Matches T == cv::Point, Point2i, etc
                raw_angle = std::atan2 (pt.y, pt.x);
//          } else if constexpr (std::something to identify morph::vVector or morph::Vector or BezCoord etc) {
            } else {
                raw_angle = std::atan2 (pt[0], pt[1]);
            }

            this->wind (raw_angle);
        }

        //! Update this->angle, this->angle_last and this->angle_sum based on \a raw_angle
        void wind (const double& raw_angle) {

            // Convert -pi -> 0 -> +pi range of atan2 to 0->2pi:
            this->angle = raw_angle >= 0 ? raw_angle : (morph::TWO_PI_D + raw_angle);

            // Set the initial angle.
            if (this->angle_last == double{-100.0}) {
                this->angle_last = angle;
            }

            double delta = double{0.0}; // delta is 'angle change'
            if (this->angle == double{0.0}) {
                // Special treatment
                if (this->angle_last > morph::PI_D) {
                    // Clockwise to 0
                    delta = (morph::TWO_PI_D - this->angle_last);
                } else if (this->angle_last < morph::PI_D) {
                    // Anti-clockwise to 0
                    delta = -this->angle_last;
                } else { //angle_last must have been 0.0
                    delta = double{0.0};
                }

            } else {

                // Special treatment required ALSO if we crossed the 0 line without being on it.
                if (this->angle_last > morph::PI_D && this->angle < morph::PI_D) {
                    // crossed from 2pi side to 0 side: Clockwise
                    delta = this->angle + (morph::TWO_PI_D - this->angle_last);
                } else if (this->angle_last < morph::PI_OVER_2_D && this->angle > morph::PI_x3_OVER_2_D) {
                    // crossed from 0 side to 2pi side: Anti-clockwise
                    delta = - this->angle_last - (morph::TWO_PI_D - this->angle);
                } else { // Both are > pi or both are < pi.
                    delta = (this->angle - this->angle_last);
                }
            }
            this->angle_last = this->angle;
            this->angle_sum += delta;
        }

        //! Reset the angle member variables
        void reset() {
            this->angle_last = double{-100.0};
            this->angle = double{0.0};
            this->angle_sum = double{0.0};
        }

        //! Member reference to the boundary
        const Container<T, Alloc>& boundary;
        //! Current angle around a point
        double angle;
        //! The sum of angles
        double angle_sum;
        //! The angle of the last boundary point
        double angle_last;
    };

} // namespace morph
