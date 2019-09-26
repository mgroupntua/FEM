using System;
using System.Collections.Generic;
using System.IO;

namespace MGroup.FEM.Readers
{
	enum NEUPosition
	{
		Undefined = 0,
		SectionStart,
		Section,
		SectionEnd
	}

	public class NEUReader
	{
		private const string sectionSeparator = "   -1";
		private readonly string fileName;

		public NEUReader(string fileName)
		{
			this.fileName = fileName;
		}

		public IEnumerable<int> GetSections()
		{
			var sections = new List<int>();
			NEUPosition position = NEUPosition.Undefined;
			using (var sr = new StreamReader(fileName))
			{
				while (!sr.EndOfStream)
				{
					var line = sr.ReadLine();
					switch (position)
					{
						case NEUPosition.Undefined:
							if (line == sectionSeparator)
								position = NEUPosition.SectionStart;
							else
								throw new InvalidDataException("File does not start with section separator.");
							break;
						case NEUPosition.SectionStart:
							int section = 0;
							if (Int32.TryParse(line, out section))
							{
								position = NEUPosition.Section;
								sections.Add(section);
							}
							else
								throw new InvalidDataException(String.Format("Expected section number but got {0}.", line));
							break;
						case NEUPosition.Section:
							if (line == sectionSeparator)
								position = NEUPosition.SectionEnd;
							break;
						case NEUPosition.SectionEnd:
							if (line == sectionSeparator)
								position = NEUPosition.SectionStart;
							else
								throw new InvalidDataException("New section does not start with section separator.");
							break;
					}
				}
			}

			return sections;
		}
	}
}
