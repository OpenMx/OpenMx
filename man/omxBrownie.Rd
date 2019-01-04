%
%   Copyright 2007-2019 by the individuals mentioned in the source code history
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{omxBrownie}
\alias{omxBrownie}

\title{Make Brownies in OpenMx}

\description{
  This function returns a brownie recipe.
}

\usage{
omxBrownie(quantity=1, walnuts=TRUE, wfpb=FALSE)
}

\arguments{
  \item{quantity}{Number of batches of brownies desired. Defaults to one.}
  \item{walnuts}{Logical. Indicates whether walnuts are to be included in the brownies. Defaults to TRUE.}
  \item{wfpb}{Logical. Indicates whether to display the whole food plant based version.  Defaults to FALSE.}
}

\details{
    Returns a brownie recipe. Alter the `quantity` variable to make more
    pans of brownies. Ingredients, equipment and procedure are listed,
    but neither ingredients nor equipment are provided.
    
    Raw cocoa powder can be used instead of Dutch processed cocoa
    for approximately double the antioxidants and flavonols. However,
    raw cocoa powder is not as smooth and delicious in taste.
    
    For the whole food plant based (wfpb) version of the recipe,
    we substitute coconut butter for dairy butter because dairy butter
    contains a large proportion of saturated fat that raises deadly LDL
    cholesterol (Trumbo & Shimakawa, 2011). In contrast, coconut butter
    has so much fiber that the considerable saturated fat that it
    contains is mostly not absorbed (Padmakumaran, Rajamohan & Kurup,
    1999).
    You can substitute erythritol (den Hartog et al, 2010) for sucanat
    (Lustig, Schmidt, & Brindis, 2012) to improve the glycemic index and
    reduce calorie density. We substitute whole wheat flour for all-purpose
    wheat flour because whole grains are associated with improvement in
    blood pressure (Tighe et al, 2010).  }

\value{
  Returns a brownie recipe. 
}

\references{

Padmakumaran Nair K.G, Rajamohan T, Kurup P.A. (1999). Coconut kernel protein
modifies the effect of coconut oil on serum lipids. Plant Foods Hum
Nutr. 53(2):133-44.

Tighe P, Duthie G, Vaughan N, Brittenden J, Simpson W.G, Duthie S, Mutch
W, Wahle K, Horgan G, Thies F. (2010). Effect of increased consumption of
whole-grain foods on blood pressure and other cardiovascular risk
markers in healthy middle-aged persons: a randomized controlled
trial. Am. J. Clin. Nutr. 92(4), 733-40.

R H Lustig, L A Schmidt, C D Brindis. (2012). Public health: The toxic truth
about sugar. Nature. 482 27-29.

den Hartog G.J, Boots A.W, Adam-Perrot A, Brouns F, Verkooijen I.W, Weseler
A.R, Haenen G.R, Bast A. (2010). Erythritol is a sweet 
antioxidant. Nutrition. 26(4), 449-58.

Trumbo P.R, Shimakawa T. (2011). Tolerable upper intake levels for trans
fat, saturated fat, and cholesterol. Nutr. Rev. 69(5), 270-8. doi:
10.1111/j.1753-4887.2011.00389.x.

The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\seealso{
More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}

\examples{
# Return a brownie recipe
omxBrownie()
}
