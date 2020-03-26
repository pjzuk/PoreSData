/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fixedFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    nName_("n"),
    snPhiGradName_("snPhiGrad"),
    sign_(1.0),
    D_(1.0),
    mu_(100.0)
{}

Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fixedFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    nName_(ptf.nName_),
    snPhiGradName_(ptf.snPhiGradName_),
    sign_(ptf.sign_),
    D_(ptf.D_),
    mu_(ptf.mu_)
{}

Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),    
    nName_(dict.lookup("n")),
    snPhiGradName_(dict.lookupOrDefault<word>("snPhiGrad", "snPhiGrad")),
    sign_(readScalar(dict.lookup("sign"))),
    D_(readScalar(dict.lookup("D"))),
    mu_(readScalar(dict.lookup("mu")))
{
    if (dict.found("gradient")) 
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}

Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fixedFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    nName_(wbppsf.nName_),
    snPhiGradName_(wbppsf.snPhiGradName_),
    sign_(wbppsf.sign_),
    D_(wbppsf.D_),
    mu_(wbppsf.mu_)
{}

Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fixedFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    nName_(wbppsf.nName_),
    snPhiGradName_(wbppsf.snPhiGradName_),
    sign_(wbppsf.sign_),
    D_(wbppsf.D_),
    mu_(wbppsf.mu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& np =
        patch().lookupPatchField<volScalarField, vector>(nName_);

    const fvsPatchField<scalar>& snPhiGrad =
        patch().lookupPatchField<surfaceScalarField, scalar>(snPhiGradName_);

    gradient() =  - sign_ * np * (mu_/ D_ ) * (snPhiGrad / patch().magSf()) ;

    fixedGradientFvPatchScalarField::updateCoeffs();

}


void Foam::fixedFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value",os);
    writeEntryIfDifferent<word>(os, "snPhiGrad", "snPhiGrad", snPhiGradName_);
    os.writeKeyword("n") << nName_ << token::END_STATEMENT << nl;
    os.writeKeyword("sign") << sign_ << token::END_STATEMENT << nl;
    os.writeKeyword("D") << D_ << token::END_STATEMENT << nl;
    os.writeKeyword("mu") << mu_ << token::END_STATEMENT << nl;
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFluxFvPatchScalarField
    );
}

// ************************************************************************* //
