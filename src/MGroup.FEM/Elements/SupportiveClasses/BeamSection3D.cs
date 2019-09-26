namespace MGroup.FEM.Elements.SupportiveClasses
{
	/**
     * Class for three-dimensional beam element sections.
     *
     * @author Theofilos Manitaras
     */
	public class BeamSection3D
	{
		private readonly double torsionalInertia;
		private readonly double area;
		private readonly double inertiaZ;
		private readonly double effectiveAreaY;
		private readonly double effectiveAreaZ;
		private readonly double inertiaY;

		/**
         * Creates a new instance of {@link BeamSection3D} class.
         *
         * @param area
         *            Section area
         * @param inertiaY
         *            Moment of inertia on axis Y
         * @param inertiaZ
         *            Moment of inertia on axis Z
         * @param torsionalInertia
         *            Torsional moment of inertia
         * @param effectiveAreaY
         *            The effective area on axis Y
         * @param effectiveAreaZ
         *            The effective area on axis Z
         */
		public BeamSection3D(double area, double inertiaY, double inertiaZ, double torsionalInertia, double effectiveAreaY, double effectiveAreaZ)
		{
			this.area = area;
			this.inertiaY = inertiaY;
			this.inertiaZ = inertiaZ;
			this.torsionalInertia = torsionalInertia;
			this.effectiveAreaY = effectiveAreaY;
			this.effectiveAreaZ = effectiveAreaZ;
		}

		public double Area { get { return area; } }
		public double EffectiveAreaY { get { return effectiveAreaY; } }
		public double EffectiveAreaZ { get { return effectiveAreaZ; } }
		public double InertiaY { get { return inertiaY; } }
		public double InertiaZ { get { return inertiaZ; } }
		public double TorsionalInertia { get { return torsionalInertia; } }
	}
}
