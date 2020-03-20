using System.Runtime.InteropServices;
using Xunit;

// General Information about an assembly is controlled through the following 
// set of attributes. Change these attribute values to modify the information
// associated with an assembly.

// Setting ComVisible to false makes the types in this assembly not visible 
// to COM components.  If you need to access a type in this assembly from 
// COM, set the ComVisible attribute to true on that type.
[assembly: ComVisible(false)]

// The following GUID is for the ID of the typelib if this project is exposed to COM
[assembly: Guid("33e1d346-f554-4db9-8017-9305b201cb64")]

// Version information for an assembly consists of the following four values:
//
//      Major Version
//      Minor Version 
//      Build Number
//      Revision
//

// The following will prevent XUnit from running tests in parallel. Normally it should not be necessary, but up until 10/12/2018
// it was not worth the effort to make the classes used in some tests thread-safe.
[assembly: CollectionBehavior(DisableTestParallelization = true)]